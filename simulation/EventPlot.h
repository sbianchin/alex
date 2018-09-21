#ifndef EVENT_PLOT_H
#define EVENT_PLOT_H

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2F.h"

class EventPlot {
	public:
		EventPlot(const char* title);

		TGraphErrors* ShowBars(std::vector<int> bars, Color_t color=kBlack);
		TGraphErrors* ShowTGraphErrors(TGraphErrors* graph);
		TF1* ShowTF1(TF1* func_plots, int color=-1);
		TH2F* ShowPoint(Double_t x, Double_t y, Color_t color=kBlack, Style_t=kMultiply, const char* name="");
		TArrow* ShowTArrow(TArrow* arrow, Color_t color=kBlack);
		void ShowValue(const char* name, Double_t value);
		void ShowValue(const char* name, Int_t value);

		void Draw(TCanvas* canvas, unsigned int subplot);

	private:
		/*
		 * Tagged value type to hold either Double_t or Int_t and stores the
		 * type information in a enum. Used to store values of both types in
		 * one map and adjust printing accordingly to type.
		 */
		struct Value {
			enum {
				Int,
				Double
			} type;

			union {
				Int_t i;
				Double_t d;
			};
		};

		const char* title;

		std::vector<TGraphErrors*> graph_plots;
		std::vector<TF1*> func_plots;
		std::vector<TH2F*> point_plots;
		std::vector<TArrow*> arrow_plots;
		std::map<const char*, Value> value_map;
};

#endif
