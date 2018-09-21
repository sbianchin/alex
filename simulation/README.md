# Single Event Analysis
A single event analysis is performed by loading the macro `Many_Tracks_Fit_2.2.C` and then calling the function `FindTwoTracks(int run_number, int event_number)`.
The first argument is MC run number and the second the event number.
The data file is expected to be named `G4Run{run_number}.root` and has to be located in the same directory.
It will open an `EventPlot` from `EventPlot.C` to show the steps of the analysis procedure.

    $ root
	root [0] .L Many_Tracks_Fit_2.2.C+
	root [1] FindTwoTracks(828, 5)

The full signature of the function reads

	Lepton FindTwoTracks(Int_t Run_Number, Int_t ievt=0, int savefile = 0, Int_t Switch_Display = 1, Int_t Switch_Output = 1, Int_t batch = 0);

If `savefile` is set to true, the result will be written in CSV format to the file `RUN_{run_number}_Find_Two_Tracks.txt`.
If `Switch_Display` is set to false, it will not produce a graphical display.
If `Switch_Output` is set to false, output to the standard output stream `cout` (`stdout`) is suppressed.
If `batch` is set to true, the graphical output will be redirected into a PDF.
This is mostly used for the batch analysis called from another program.

Finding the kaon and lepton tracks is done in `FindKaons.C`, the two track separation in `Many_Tracks_Fit_2.2.C` in the function `FindTrack()`.
Parameters of the Geant4 simulation are set it `G4DataRootApr19.h`, general parameters in `GlobalParameters.h` and those to be varied systematically in `TweakParameters.h`.

# Batch Analysis
## Analysing Whole Data Sets
To analyse multiple events in one MC run file, load the macro `Batch_Two_Tracks_Fit.C` and call the function `Batch_Two_Tracks_Fit(int run_number, int flag)`.
The data file is expected to be named `G4Run{run_number}.root` and has to be located in the same directory.
This will analyze the whole data set and save a CSV file of output variables including a header to `RUN_{run_number}_Find_Two_Tracks.txt`.
To analyse only parts of the data set, set the third and fourth parameter to the minimum and maximum event number, respectively.
If the flag is set to false, only text output is produced, if it is set to false, only text output is produced.
If it is set to true and the number of events is less than 1000, the graphical output of `FindTwoTracks()` is saved in a PDF.

    $ root
	root [0] .L Batch_Two_Tracks_Fit.C+
	root [1] Batch_Two_Tracks_Fit(828,0) // Generates CSV for whole data file
	root [2] Batch_Two_Tracks_Fit(828,1,50,100) // Saves plots for events 50 to 100

It can also be called directly from the commandline

    $ root Batch_Two_Tracks_Fit.C+(828,0)

To conveniently run this on multiple data sets at once, call the function `Multi_Batch(int count, ...)`.
The first argument is the total number of data files to be analyzed, the rest is an arbitrary number of run numbers.

## Meta-analysis
To run the meta analysis Python program, you need to install Python 3 and the following libraries via `pip`:

    $ pip3 install scipy pandas matplotlib

You can either load the file `paramter_study.py` into iPython or Jupyter to use the functions interactively or edit the file and call them directly.


`analyze_run(run_num, is_filename=False)` is used to analyze a single data file from `Batch_Two_Tracks_Fit()`. 
The data file has to be in the same directory and is called `RUN_{run_num}_Find_Two_Tracks.txt`.
If `is_filename` is set to `True`, an arbitrary file name can be given. 
This is designed to be used in interactive mode, like in jupyter or iPython.
If the lines 68 to 71 are made active again by removing the `#`, a file `problematic_events.txt` will be procuded which lists events with bad primary angle error and K-stop.

`plot_analysis(runs, filename)` can be used to run `analyze_run()` on the array `runs` of run_numbers and save all summary plots in the PDF `filename`.

	[python file parameter_study.py]
	plot_analysis(828, "two-track-batch-meta_{}.pdf".format(runs))

	[commandline]
	$ python3 paramter_study.ph

# Parameter Study
In order to vary parameters, `parameter_study(run_numbers, tweak_parameters_variations, from_scratch=False, run_labels=None)`.
It takes an array of run numbers and a dictionary of mappings of parameter names to arrays or ranges of values.
This will generate a `TweakParameters.h` file with the corresponding set of parameter values for each variation, data file and parameter.
The root will be called to run the batch analysis and the results are combined in CSV and PDF files, sorted by parameter and run file.
The default values that should be used for all the other parameters apart from the one that is varied, are set in the function, from line 263 on.
If `from_scratch` is set to `True`, the root analysis is carried out for every variation, even if the corresponding CSV exists already.
Otherwise, the cached data will be used.
To label the runs not by number but name in the legend of the final overview plot, pass a dictionary mapping run numbers to labels as `run_labels`.

	[python file parameter_study.py]
	runs = [828, 2141, 2161, 3212]
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

	parameter_study(runs, tweak_parameters_variations, run_labels=run_labels)


	[commandline]
	$ python3 paramter_study.py

This would run the parameter variation of all parameters over different ranges for the run files 828, 214-1, 216-1, and 321-2.
The will receive custom labels.
