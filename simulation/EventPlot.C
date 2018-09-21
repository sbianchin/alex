/*
 * EventPlot
 * 
 * Used to print a display a single event with different types of overlays, such
 * as hit bars, fits or point markers. It is designed to be used with a
 * predevided canvas.
 *
 * Example usage:
 *
 *   TCanvas* canvas = new TCanvas("canvas", "EventPlot Example" ,0 ,200, 1050, 700);
 *   canvas->Divide(3,2);
 *   EventPlot* event_plot("Example Plot");
 *   event_plot->ShowBars(lepton_bars, kRed);
 *   event_plot->ShowTGraphErrors(lepton_fit);
 *   event_plot->ShowTF1(kaon_fit);
 *   event_plot->ShowPoint(x, y, kBlue, kStar);
 *   event_plot->ShowTArrow(arrow);
 *   event_plot->ShowValue("Angle", angle);
 *   event_plot->Draw(canvas, 2);
 *
 * where
 *   lepton_bars  std::vector<int> of target bar indices
 *   lepton_fit   TGraphErrors object (pointer) to be drawn
 *   kaon_fit     TF1 object (pointer) function to be drawn
 *   arrow        TArrow object (pointer) to be drawn
 *   angle        real value to be printed
 *   x, y         real coordinates for a marker
 * and the plot is drawn to the 2nd subplot. This needs an Int_t tof1N[5]
 * array of hit TOF1 counters defined as a global variable in G4DataRootApr19.h
 *
 * When the Show…() methods are called, the respective overlay is added to an
 * internal list. On a Draw() call, all background objects for the general plot
 * layout are drawn, then the list of overlays is iterated and drawn. The
 * default styling of the overlays is done in the Show…() methods. Most of the
 * Show…() methods return a pointer to the object to be drawn, so that line
 * attributes can be further manipulated if one is not happy with the default
 * styling.
 *
 * The class is designed to be used like this: First create an EventPlot object
 * for the kind of display you want to generate from calculated data, possibly
 * to show intermediate steps in an algorithm or single results. Add bars, fits
 * and values into the plot along their calculation in the algorithm. In the
 * main program, create a canvas and draw all the EventPlot objects in the end.
 */

/* TODO
 * move to ShowTOF1(std::vector<int> tof1_hit) function to highlight hit TOF1s
 * add templated function that can draw anything with a Draw() method
 */


#include "EventPlot.h"

#include "TPaveText.h"
#include "TGraph.h"
#include "TEllipse.h"
#include "TGaxis.h"

#include <sstream>

using namespace std; // parameter header files need this
#include "G4DataRootApr19.h"
#include "GlobalParameters.h"
#include "TargetParameters.h"

/* Constructor for EventPlot object. 
 *
 * const char* title_  title string shown above the plot
 */
EventPlot::EventPlot(const char* title_) : title(title_) {}

/*
 * Show target bars. Takes a std::vector<int> of bar indices starting with 0
 * at the top row left and ending with 255 in the bottom row right. A bar color
 * can be passed as an optional argument. Defaults to black
 *
 * std:vector<int> bars  bar indices to be shown
 * Color_t color         optional ROOT color for the bars
 */
TGraphErrors* EventPlot::ShowBars(std::vector<int> bars, Color_t color) {
	std::vector<Double_t> coords_x;
	std::vector<Double_t> coords_y;
	std::vector<Double_t> errors_x;
	std::vector<Double_t> errors_y;

	for (int bar : bars) {
		coords_x.push_back(Xloc[bar]);
		coords_y.push_back(Yloc[bar]);
		errors_x.push_back(TARGET_Errors_X);
		errors_y.push_back(TARGET_Errors_Y);
	}

	TGraphErrors *bar_plot = new TGraphErrors(
		bars.size(),
		&coords_x[0],
		&coords_y[0],
		&errors_x[0],
		&errors_y[0]
	);

	bar_plot->GetXaxis()->SetLimits(-50.,50.);
	bar_plot->GetYaxis()->SetRangeUser(-50.,50.);
	bar_plot->SetMarkerStyle(kFullSquare);
	bar_plot->SetMarkerColor(color);
	bar_plot->SetMarkerSize(0.8);

	graph_plots.push_back(bar_plot);
	return graph_plots.back();
}

/*
 * Show a TGraphErrors object. This is to be expect to have all styles and
 * colors already set, so no draw options can be passed here. This can be used
 * to show bars and fits. Returns a pointer to the TGraphErrors object ot be
 * drawn. This can be used to manipulate the drawing attributes.
 *
 * TGraphErrors* graph  the graph to be drawn
 */
TGraphErrors* EventPlot::ShowTGraphErrors(TGraphErrors* graph) {
    graph->GetXaxis()->SetLimits(-50.,50.);
    graph->GetYaxis()->SetRangeUser(-50.,50.);
	graph_plots.push_back(graph);
	return graph_plots.back();
}

/*
 * Show a TF1 function object. Line color can be passed as optional argument.
 * Defaults to black. Returns a pointer to the TF1 object that is to be plotted.
 * This can be used to manipulate style attributes.
 *
 * TF1* func  function object to be drawn
 * int color  optional ROOT color
 */
TF1* EventPlot::ShowTF1(TF1* func, int color) {
    func->GetXaxis()->SetLimits(-50.,50.);
    func->GetYaxis()->SetRangeUser(-50.,50.);
	func->SetLineWidth(2);
	// negative default value is used to indicate it was not set at call and
	// the color of the object should be used instead
	if (-1 != color) {
		func->SetLineColor(color);
	}

	func_plots.push_back(func);
	return func_plots.back();
}

/*
 * Show a marker for a single point on the plot. Color and style can be passed
 * as optional arguments. Defaults to a black x. Returns a pointer to the TH2F
 * object of the point. This can be used to manipulate the marker styles.
 *
 * Double_t x        x coordinate, target center is at (0,0)
 * Double_t y        y coordinate, target center is at (0,0)
 * Color_t color     optional ROOT color
 * Style_t style     optional maker style from https://root.cern.ch/doc/master/classTAttMarker.html
 * const char* name  optional name of the marker, is set as a LineAttribute
 */
TH2F* EventPlot::ShowPoint(Double_t x, Double_t y, Color_t color, Style_t style, const char* name) {
	TH2F* point_plot = new TH2F(name, name, 500, -50, 50, 500, -50, 50);
	point_plot->SetStats(false);
	point_plot->Fill(x,y);
	point_plot->SetMarkerStyle(style);
	point_plot->SetMarkerColor(color);
	point_plot->SetMarkerSize(2.5);

	point_plots.push_back(point_plot);
	return point_plots.back();
}

/*
 * Show a TArrow. Color can be passed as an optional argument, black is the
 * default value. A pointer to the TArrow object to be drawn is returned for
 * style manipulation.
 *
 * TArrow* arrow  Arrow to be drawn
 * Color_t color  optional ROOT color for the arrow
 */
TArrow* EventPlot::ShowTArrow(TArrow* arrow, Color_t color) {
	arrow->SetLineWidth(2);
	arrow->SetLineColor(color);
	arrow_plots.push_back(arrow);
	return arrow_plots.back();
}

/*
 * Print a key-value-pair of a string and a numerical value as a legend beside
 * the plot. Will be printed in "name = 23" format. Currently real values as
 * Double_t and integers as Int_t are supported. Real values will be printed
 * with fixed two digits after decimal precision.
 *
 * An internal data type Value is used to store the values of both types in
 * one map, together with type information used for printing.
 * 
 * const char* name          name of the value to be printed
 * typename value_to_show    value to be printed
 */
void EventPlot::ShowValue(const char* name, Double_t value_to_show) {
	Value value;
	value.type = Value::Double;
	value.d = value_to_show;
	value_map[name] = value;
}

void EventPlot::ShowValue(const char* name, Int_t value_to_show) {
	Value value;
	value.type = Value::Int;
	value.i = value_to_show;
	value_map[name] = value;
}


/*
 * Draw the whole plot including overlays.
 *
 * TCanvas* canvas       canvas to be drawn on
 * unsigned int subplot  subplot number in the canvas to which to draw to
 */
void EventPlot::Draw(TCanvas* canvas, unsigned int subplot) {
	canvas->cd(subplot)->Range(-50, -50, 50, 50);

	// setup a coordinate grid an axes with ticks for all other plots
	TH2F *empty = new TH2F(title, title, 500, -50, 50, 500, -50, 50);
	empty->SetStats(false);
	empty->Draw();

	TGaxis *A1 = new TGaxis(-50,50,50,50,"pol1",510,"-U");
	TGaxis *A2 = new TGaxis(50,-50,50,50,"pol1",510,"+U");
	A1->Draw();
	A2->Draw();

	TH2F *target_center = new TH2F("target_center", "Target Center", 500, -50, 50, 500, -50, 50);
	target_center->SetStats(false);
	target_center->Fill(0., 0.);
	target_center->SetMarkerStyle(5);
	target_center->SetMarkerColor(1);
	target_center->SetMarkerSize(2);
	target_center->Draw("same");

	// Draw TOF1 gaps and highlight hit ones
	TLine* gaps[12] = {Gap1l, Gap2l, Gap3l, Gap4l, Gap5l, Gap6l, Gap7l, Gap8l, Gap9l,Gap10l, Gap11l, Gap12l};
	for (TLine* gap : gaps) {
		gap->SetLineWidth(10);
		gap->SetLineColor(15);
	}
	for (int i = 0; i < 5; i++) {
		if (tof1N[i] != -999) {
			gaps[tof1N[i]-1]->SetLineColor(42);
		}
	}
	for (TLine* gap : gaps) {
		gap->Draw();
	}

	const float R_TOF1 = 47.1;
	TEllipse *bound_tof1 = new TEllipse(0, 0, R_TOF1, 0);
	bound_tof1->SetFillStyle(0);
	bound_tof1->SetLineColor(kMagenta);
	bound_tof1->SetLineWidth(1);
	bound_tof1->Draw();
	const float R_TARGET = 29.0;
	TEllipse *bound_target = new TEllipse(0, 0, R_TARGET, 0);
	bound_target->SetFillStyle(0);
	bound_target->SetLineColor(kBlack);
	bound_target->SetLineWidth(1);
	bound_target->Draw();
	const float R_SFT_L1 = 40.0;
	TEllipse *bound_sft = new TEllipse(0, 0, R_SFT_L1, 0);
	bound_sft->SetFillStyle(0);
	bound_sft->SetLineColor(kBlue);
	bound_sft->SetLineWidth(1);
	bound_sft->Draw();

	for (TGraphErrors* graph_plot : graph_plots) {
		graph_plot->Draw("P");
	}

	for (TF1* func_plot : func_plots) {
		func_plot->Draw("same");
	}

	for (TH2F* point_plot : point_plots) {
		point_plot->Draw("same");
	}

	for (TArrow* arrow_plot : arrow_plots) {
		arrow_plot->Draw("same");
	}

	const int num_values = value_map.size();
	if (num_values > 0) {
		// Restrict height of the box so that it won't cover the important
		// part of the plot, but only do so above a max number of values.
		// For 4 or less values, scale the box accordingly.
		float rel_height = 0.2;
		if (num_values < 5) {
			rel_height = 0.05 * num_values;
		}
		TPaveText* text_box = new TPaveText(0.6, 1-rel_height, 1, 1, "NDC");
		for (const auto& kv : value_map) {
			std::stringstream line;
			switch (kv.second.type) {
				case Value::Double:
					line.setf(ios::fixed);
					line.precision(2);
					line << kv.first << " = " << kv.second.d;
					break;
				case Value::Int:
					line << kv.first << " = " << kv.second.i;
					break;
			}
			text_box->AddText(line.str().c_str());
		}
		text_box->Draw("same");
	}
}

