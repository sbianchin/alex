# parameter_study.py
#
# example usage is given below the function definitions.
#
# analyze_run(run_num, is_filename=False)
#   analyze a single data file. Data file has to be in the same directory and
#   is called "RUN_{run_num}_Find_Two_Tracks.txt". If `is_filename` is set to
#   True, an arbitrary file name can be given. This is designed to be used in
#   interactive mode, like in jupyter or iPython
# plot_analysis(runs, filename)
#   run `analyze_run()` on the array `runs` of run_numbers and save all
#   summary plots in the PDF `filename`.
# parameter_study(run_numbers, tweak_parameters_variations, from_scratch=False, run_labels=None)
#   vary the parameters in `TweakParameters.h` according to the dict
#   `tweak_parameters_variations` for runs in `run_numbers`. If the runs should
#   have named labels in the legend, pass the dict `run_labels` with a mapping.
#   from_scratch will cause the recalculation for all data sets, otherwise, old
#   CSV files are reused. Default values are defined in the function definition.
#
# projection_plot(ax, frame, column_x, column_y, frame_one_track=None, commands=None, bins=100)
#   helper function to produce scatter plots with histogram projections for
#   both axis.

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams['legend.framealpha'] = None
mpl.rcParams['legend.fancybox'] = False
mpl.rcParams['legend.numpoints'] = 1
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cycler import cycler

import os
import shutil
from collections import defaultdict
import subprocess


def analyze_run(run_num, is_filename=False):
    if is_filename:
        data_file = run_num
    else:
        data_file = 'RUN_{}_Find_Two_Tracks.txt'.format(run_num)
        if not os.path.isfile(data_file):
            cmd = ['root', '-e', '.L Batch_Two_Tracks.C', '-e', 'Batch_Two_Tracks_Fit({}, 0);'.format(run_num)]
            #subprocess.Popen(cmd).wait()

    data = np.genfromtxt(data_file, delimiter=',')
    frame = pd.DataFrame(data=data, index=data[:,1].astype(int), columns=["Run_Number", "ievt", "reduced_ChiS", "ndf", "first_track_ChiS_reduced", "ndf_2", "second_track_ChiS_reduced", "ndf_3", "E_positron", "delta_energy", "delta_length", "delta_length_xy", "angle_between", "two_track_fit", "angle_primary_sim", "angle_lepton_all" ,"angle_primary", "delta_angle_primary", "angle_secondary", "delta_angle_secondary", "k_stop_error", "angle_primary_error", "targL", "no_leptons", "no_kaons", "k_stop_radius", "TWO_TRACK_MIN_CHISQ", "TWO_TRACK_MIN_LEPTON_BARS", "FIT_TOF1_WEIGHT", "K_STOP_CENTROID_THRESH", "PATH_TRAVERSIAL_USE_ALL", "PATH_TRAVERSIAL_DIJKSTRA_JUMP_RADIUS", "PATH_TRAVERSIAL_ALL_PENALTY"])

    cuts = []
    cut_no_empty_tracks              = ( (frame['no_leptons'] == 0) & (frame['no_kaons'] == 0) )
    frame = frame[cut_no_empty_tracks]
    cut_k_stop_error_gt_10           = (frame['k_stop_error']        > 10)
    cut_angle_primary_error_positive = (frame['angle_primary_error'] > 0)
    cuts.append(cut_angle_primary_error_positive)
    cut_k_stop_error_lt_100          = (frame['k_stop_error']        < 100)
    cuts.append(cut_k_stop_error_lt_100)
    cut_angle_primary_error_lt_10    = (frame['angle_primary_error'] < 10)
    cuts.append(cut_angle_primary_error_lt_10)
    cut_faulty_delta_energy          = (frame['delta_energy'] > 0)
    cuts.append(cut_faulty_delta_energy)
    cut_large_delta_energy           = (frame['delta_energy'] < 120)
    cuts.append(cut_large_delta_energy)

    #with open("problematic-events.txt", 'w') as file:
    #    file.write(frame[cut_k_stop_error_gt_10].sort_values(by='k_stop_error', ascending=False)['k_stop_error'].to_csv(sep="\t", header=True))
    #    file.write("\n\n")
    #    file.write(frame[~cut_angle_primary_error_lt_10].sort_values(by='angle_primary_error', ascending=False)['angle_primary_error'].to_csv(sep="\t", header=True))

    for cut in cuts:
        frame = frame[cut]

    num_events = len(frame.index)
    delta_detection_efficiency = np.count_nonzero(frame['two_track_fit'] == 1) / num_events


    results = [
            int(frame['Run_Number'].mean()),
            len(frame['ievt']),
            delta_detection_efficiency,
            frame['angle_primary_error'].mean(),
            frame['angle_primary_error'].std(),
            np.nanmean(frame['delta_angle_primary']),
            frame['k_stop_error'].mean(),
            cut_k_stop_error_gt_10.sum()/num_events,
            (~cut_angle_primary_error_lt_10).sum()/num_events,
            (~cut_no_empty_tracks).sum()/num_events,
            ]

    info = []
    info.append(["run number", "{}".format(int(frame['Run_Number'].mean()))])
    info.append(["number of events", "{}".format(len(frame['ievt']))])
    info.append(["delta detection efficiency", "{:.2f}%".format(delta_detection_efficiency*100)])
    info.append(["mean primary angle error", "{:.2f}°".format(frame['angle_primary_error'].mean())])
    info.append(["primary angle error std deviation", "{:.2f}°".format(frame['angle_primary_error'].std())])
    info.append(["mean difference one/two-track fit", "{:.2f}°".format(np.nanmean(frame['delta_angle_primary']))])
    info.append(["mean K-stop error", "{:.2f}mm".format(frame['k_stop_error'].mean())])
    info.append(["fraction events k_stop_error > 10", "{:.2f}%".format(cut_k_stop_error_gt_10.sum()/num_events*100)])
    info.append(["fraction events angle_primary_error > 10", "{:.2f}%".format((~cut_angle_primary_error_lt_10).sum()/num_events*100)])
    info.append(["events with missing tracks (kaon/lepton)", "{:.2f}%".format((~cut_no_empty_tracks).sum()/num_events*100)])

    cut_two_track_fit = (frame['two_track_fit'] == 1)
    one_track_frame = frame[~cut_two_track_fit]
    frame = frame[cut_two_track_fit]
    #print(frame.sort_values(by='angle_primary_error', ascending=False)[['two_track_fit','angle_primary_error']].head())
    #print(frame.sort_values(by='delta_angle_primary', ascending=False)[['two_track_fit','delta_angle_primary']].head())
    #print(frame.sort_values(by='delta_energy', ascending=False)[['angle_primary_error', 'targL', 'delta_energy']].head(10))

    PLOT_ROWS=3
    PLOT_COLS=4
    PLOT_SIZE=6
    fig = plt.figure(figsize=(PLOT_COLS*PLOT_SIZE, PLOT_ROWS*PLOT_SIZE))
    plt.subplots_adjust(wspace=0.3, hspace=0.3)

    plt.subplot(PLOT_ROWS, PLOT_COLS, 2)
    plt.axis('off')
    info = np.array(info)
    info_table = plt.table(rowLabels=info[:,0], cellText=info[:,1].reshape(-1,1), colWidths=[0.2], loc='center left')
    info_table.set_fontsize(14)
    info_table.scale(2,1.5)

    plt.subplot(PLOT_ROWS, PLOT_COLS, 4)
    plt.axis('off')
    params = ["TWO_TRACK_MIN_CHISQ", "TWO_TRACK_MIN_LEPTON_BARS", "FIT_TOF1_WEIGHT", "K_STOP_CENTROID_THRESH", "PATH_TRAVERSIAL_USE_ALL", "PATH_TRAVERSIAL_DIJKSTRA_JUMP_RADIUS", "PATH_TRAVERSIAL_ALL_PENALTY"]
    params_t = []
    for param in params:
        params_t.append([param, "{}".format(frame[param].mean())])
    params_t = np.array(params_t)
    params_table = plt.table(rowLabels=params_t[:,0], cellText=params_t[:,1].reshape(-1,1), colWidths=[0.2], loc='center left')
    params_table.set_fontsize(14)
    params_table.scale(2,1.5)

    projection_plot(
        plt.subplot(PLOT_ROWS, PLOT_COLS, PLOT_COLS+1),
        frame, "angle_between", "angle_primary_error", commands=lambda ax: (
        ax.legend(),
        ax.set_yscale("log"),
        ax.set_title("Angle error by angle between primary and delta", y=1.4),
        ax.set_xlabel("angle_between in degrees"),
        ax.set_ylabel("angle_primary_error in degrees")))

    projection_plot(
        plt.subplot(PLOT_ROWS, PLOT_COLS, PLOT_COLS+2),
        frame, "targL", "angle_primary_error", one_track_frame, commands=lambda ax: (
        ax.legend(),
        ax.set_title("Angle error by primary track length in target", y=1.4),
        ax.set_xlabel("track length in mm"),
        ax.set_ylabel("angle_primary_error in degrees")))

    projection_plot(
        plt.subplot(PLOT_ROWS, PLOT_COLS, PLOT_COLS+3),
        frame, "k_stop_error", "angle_primary_error", one_track_frame, commands=lambda ax: (
        ax.legend(),
        ax.set_xscale('log'),
        ax.set_xlim(left=10e-3),
        ax.set_yscale('log'),
        ax.set_ylim(bottom=10e-3),
        ax.set_title("Angle error by K-stop error", y=1.4),
        ax.set_xlabel("k_stop_error in mm"),
        ax.set_ylabel("angle_primary_error in degrees")))

    projection_plot(
        plt.subplot(PLOT_ROWS, PLOT_COLS, PLOT_COLS+4),
        frame, "delta_energy", "angle_primary_error", commands=lambda ax: (
        ax.legend(),
        ax.set_title("Angle error by delta energy", y=1.4),
        ax.set_xlabel("delta_energy in MeV"),
        ax.set_ylabel("angle_primary_error in degrees")))

    plt.subplot(PLOT_ROWS, PLOT_COLS, 2*PLOT_COLS+1)
    plt.title("Angle error distribution")
    plt.xlabel("angle_primary_error in degrees")
    plt.ylabel("number of events")
    plt.hist(frame['angle_primary_error'], bins=100, label="two track")
    plt.hist(one_track_frame['angle_primary_error'], bins=100, label="one track")
    plt.legend()
    # maybe go to two axis with normalized scaling?

    plt.subplot(PLOT_ROWS, PLOT_COLS, 2*PLOT_COLS+2)
    plt.title("One/two track fit angle difference distribution")
    plt.xlabel("delta_angle_primary in degrees")
    plt.ylabel("number of events")
    plt.hist(frame['delta_angle_primary'], bins=60, range=(frame['delta_angle_primary'].min(), 20), label="two tracks")
    plt.legend()

    plt.subplot(PLOT_ROWS, PLOT_COLS, 2*PLOT_COLS+3)
    plt.title("K-stop error distribution")
    plt.xlabel("k_stop_error in mm")
    plt.ylabel("number of events")
    plt.hist(frame['k_stop_error'], bins=100, label="two track")
    plt.hist(one_track_frame['k_stop_error'], bins=100, label="one track")
    plt.legend()

    plt.subplot(PLOT_ROWS, PLOT_COLS, 2*PLOT_COLS+4)
    plt.title("K-stop radius distribution")
    plt.xlabel("k_stop_radius in mm")
    plt.ylabel("number of events")
    plt.hist(frame['k_stop_radius'], bins=30, label="two track")
    plt.hist(one_track_frame['k_stop_radius'], bins=50, label="one track")
    plt.legend()

    plt.tight_layout()
    return fig, results

def projection_plot(ax, frame, column_x, column_y, frame_one_track=None, commands=None, bins=100):
    #ax.tick_params(top=True, left=True, bottom=True, right=True, labeltop=True, labelleft=True, labelbottom=True, labelright=True)
    ax.scatter(frame[column_x], frame[column_y], label="one track", marker='x')

    divider = make_axes_locatable(ax)
    ax_hist_x = divider.append_axes("top", 1.2, pad=0.1, sharex=ax)
    for tl in ax_hist_x.get_xticklabels():
        tl.set_visible(False)
    ax_hist_x.hist(frame[column_x], bins=50)

    ax_hist_y = divider.append_axes("right", 1.2, pad=0.1, sharey=ax)
    for tl in ax_hist_y.get_yticklabels():
        tl.set_visible(False)
    ax_hist_y.hist(frame[column_y], orientation='horizontal', bins=bins)

    if frame_one_track is not None:
        ax.scatter(frame_one_track[column_x], frame_one_track[column_y], label="two tracks", marker='.')
        ax_hist_x.hist(frame_one_track[column_x], bins=bins)
        ax_hist_y.hist(frame_one_track[column_y], orientation='horizontal', bins=bins)

    if commands is not None:
        commands(ax)

def plot_analysis(runs, filename):
    with PdfPages(filename) as pdf:
        for run in runs:
            fig, results = analyze_run(run)
            pdf.savefig(fig)
        pdf.infodict()['Title'] = "Two Track Batch Meta Analysis {}".format(runs)


def parameter_study(run_numbers, tweak_parameters_variations, from_scratch=False, run_labels=None):
    root_dir = './'
    tweak_parameters_header_file = root_dir+'TweakParameters.h'
    batch_data_file_name = root_dir+'RUN_{}_Find_Two_Tracks.txt'

    tweak_parameters_header = """
                            #ifndef TWEAK_PARAMETERS_H
                            #define TWEAK_PARAMETERS_H

                            const Double_t TWO_TRACK_MIN_CHISQ = {TWO_TRACK_MIN_CHISQ};
                            const Double_t TWO_TRACK_MIN_LEPTON_BARS = {TWO_TRACK_MIN_LEPTON_BARS};

                            const Double_t FIT_TOF1_WEIGHT = {FIT_TOF1_WEIGHT};

                            const Double_t K_STOP_CENTROID_THRESH = {K_STOP_CENTROID_THRESH}; // mm

                            const Bool_t   PATH_TRAVERSIAL_USE_ALL = {PATH_TRAVERSIAL_USE_ALL}; // use FindShortestPath() instead of dijkstra()
                            const Double_t PATH_TRAVERSIAL_DIJKSTRA_JUMP_RADIUS = {PATH_TRAVERSIAL_DIJKSTRA_JUMP_RADIUS}; // mm
                            const Double_t PATH_TRAVERSIAL_ALL_PENALTY = {PATH_TRAVERSIAL_ALL_PENALTY};

                            #endif
                            """

    #""" Default values before optimization
    tweak_parameters_default = {
            'TWO_TRACK_MIN_CHISQ': 1,
            'TWO_TRACK_MIN_LEPTON_BARS': 5,
            'FIT_TOF1_WEIGHT': 3,
            'K_STOP_CENTROID_THRESH': 8,
            'PATH_TRAVERSIAL_USE_ALL': 'true',
            'PATH_TRAVERSIAL_DIJKSTRA_JUMP_RADIUS': 9,
            'PATH_TRAVERSIAL_ALL_PENALTY': 1000
            }
    """ Optimized values
    tweak_parameters_default = {
            'TWO_TRACK_MIN_CHISQ': 1.5,
            'TWO_TRACK_MIN_LEPTON_BARS': 2,
            'FIT_TOF1_WEIGHT': 5,
            'K_STOP_CENTROID_THRESH': 2,
            'PATH_TRAVERSIAL_USE_ALL': 'false',
            'PATH_TRAVERSIAL_DIJKSTRA_JUMP_RADIUS': 7.5,
            'PATH_TRAVERSIAL_ALL_PENALTY': 1000
            }
    """

    shutil.move(tweak_parameters_header_file, tweak_parameters_header_file+'.default')

    for param, variation in tweak_parameters_variations.items():
        print("\n", param, variation)

        if not os.path.exists(param):
            os.makedirs(param)
        elif from_scratch:
            shutil.rmtree(param)
            os.makedirs(param)

        figures = defaultdict(list)
        results = []

        for value in variation:
            print('  {}'.format(value))
            tweak_parameters = tweak_parameters_default.copy()
            tweak_parameters[param] = value

            if param.startswith('PATH_TRAVERSIAL_DIJKSTRA_'):
                tweak_parameters['PATH_TRAVERSIAL_USE_ALL'] = 'false'
            elif param.startswith('PATH_TRAVERSIAL_ALL_'):
                tweak_parameters['PATH_TRAVERSIAL_USE_ALL'] = 'true'

            with open(tweak_parameters_header_file, 'w') as header_file:
                header_file.write(tweak_parameters_header.format(**tweak_parameters))

            print('    ', end='')
            for run_num in run_numbers:
                print(run_num, end=' ', flush=True)

                data_file = '{}/{}.{}={}.csv'.format(param,run_num,param,value)

                if not os.path.exists(data_file):
                    with open(os.devnull) as sink:
                        root_command = ['root', "Batch_Two_Tracks_Fit.C+({},0)".format(run_num), '-l', '-q']
                        print("\n      running `{}`…\n    ".format(' '.join(root_command)), end='')
                        subprocess.call(root_command, cwd=root_dir, stdout=sink, stderr=sink)

                    shutil.copyfile(batch_data_file_name.format(run_num), data_file)

                fig, result = analyze_run(data_file, is_filename=True)
                result.append(value)

                results.append(result)
                figures[run_num].append(fig)
            print()

        for run_num, figs in figures.items():
            with PdfPages('{}/{}.{}.pdf'.format(param,run_num,param), 'a') as pdf:
                pdf.infodict()['Title'] = "Parameter Study {} [{}]".format(param, run_num)
                for fig in figs:
                    pdf.savefig(fig)
                    plt.close(fig)

        np.savetxt('{}/{}.csv'.format(param,param), results, delimiter=",")

        columns=['Run_Number', 'num_events', 'delta_detection_efficiency', 'mean_angle_primary_error', 'angle_primary_error_std_deviation', 'mean_delta_angle_primary', 'mean_k_stop_error', 'k_stop_error_gt_10', 'angle_primary_error_gt_10', 'empty_tracks']
        columns.append(param)
        frame = pd.DataFrame(data=results, columns=columns)

        def plot_variation(subplot, title, index):
            ax = plt.subplot(subplot)
            colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
            ax.set_prop_cycle(cycler('marker', ['x', 'o', 'v', '^', '*', '+']) +
                              cycler('color', colors[:6])
                             )
            plt.title(title)
            plt.xlabel(param)
            plt.ylabel(index)
            for run_num in run_numbers:
                run_frame = frame[frame['Run_Number'] == run_num]
                if run_labels is not None and run_num in run_labels:
                    label = run_labels[run_num]
                else:
                    label = 'run {}'.format(run_num)
                plt.plot(run_frame[param], run_frame[index], linestyle='-', label=label)
            plt.legend()

        plt.figure(figsize=(13,8.5))
        plt.suptitle("{} variation".format(param))

        plot_variation(231, "Mean Primary Angle Error", "mean_angle_primary_error")
        plot_variation(232, "Primary Angle Error Standard Deviation", "angle_primary_error_std_deviation")
        plot_variation(233, "Fraction of Events w/ Primary Angle Error > 10°", "angle_primary_error_gt_10")
        plot_variation(236, "Delta Detection Efficiency", "delta_detection_efficiency")
        plot_variation(234, "Mean K-Stop Error", "mean_k_stop_error")
        plot_variation(235, "Fraction of Events w/ K-Stop Error > 10mm", "k_stop_error_gt_10")

        plt.tight_layout()
        plt.subplots_adjust(top=0.9)
        plt.savefig("{}/{}.pdf".format(param,param))

    shutil.move(tweak_parameters_header_file+'.default', tweak_parameters_header_file)


runs = [828]
run_labels = {
        828:  'only deltas',
        2141: 'mixed',
        2161: 'only deltas',
        3212: 'no deltas',
        }

tweak_parameters_variations = {
        'TWO_TRACK_MIN_CHISQ': np.arange(0.5,3,0.5),
        'TWO_TRACK_MIN_LEPTON_BARS': np.arange(0,7,1),
        'FIT_TOF1_WEIGHT': np.arange(0,7,1),
        'K_STOP_CENTROID_THRESH': np.arange(0,5,0.5),
        'PATH_TRAVERSIAL_DIJKSTRA_JUMP_RADIUS': np.arange(3,21,1.5),
        'PATH_TRAVERSIAL_ALL_PENALTY': [0, 3, 5, 10, 50],
        }

#plot_analysis(runs, "two-track-batch-meta_{}.pdf".format(runs))

parameter_study(runs, tweak_parameters_variations, run_labels=run_labels)

# optimal parameters
# just run over the two different gap methods
#parameter_study(runs, { 'PATH_TRAVERSIAL_DIJKSTRA_JUMP_RADIUS': [7.5], 'PATH_TRAVERSIAL_ALL_PENALTY': [1000] }, run_labels=run_labels)

