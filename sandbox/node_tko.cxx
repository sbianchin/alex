/*
 *
 *
 */

#include <iostream>
#include <iomanip>
#include <vector>


#include "event.h"
#include "decode_tko.cxx"


DTko tko1;
DTko tko2;
DTko tko3;

int tko1_tag_event;
int tko1_tag_spill;
int tko2_tag_event;
int tko2_tag_spill;
int tko3_tag_event;
int tko3_tag_spill;

int node_tko(unsigned int *node_top)
{
	struct event_header *header
		= reinterpret_cast<struct event_header *>(node_top);
	unsigned int *data
		= node_top
		+ sizeof(struct event_header)/sizeof(unsigned int);

	tko1_tag_event = -1;
	tko1_tag_spill = -1;
	tko2_tag_event = -1;
	tko2_tag_spill = -1;
	tko3_tag_event = -1;
	tko3_tag_spill = -1;

	if (header->size < (sizeof(struct event_header)/sizeof(unsigned int) + (3 * 2))) {
		std::cerr << "#E TKO Broken data : " << std::hex;
		for (int j = 0 ; j < 8 ; j++) {
			std::cerr << std::hex << *(node_top + j) << " ";
		}
		std::cerr << std::dec << std::endl;
	} else {
		int ev_num1;
		int ev_tag1;
		int ev_num2;
		int ev_tag2;
		int ev_num3;
		int ev_tag3;

		tko1.set_data(data + 1);
		ev_num1 = tko1.get_event_number(2);
		tko1_tag_event = ev_tag1 = tko1.get_event_tag(2);
		tko1_tag_spill = tko1.get_spill_tag(2);

		/*
		std::cout << "#D " << std::hex
			<< *(datapoint + tko.get_nword() + 0) << " "
			<< *(datapoint + tko.get_nword() + 1) << " "
			<< *(datapoint + tko.get_nword() + 2) << " "
			<< *(datapoint + tko.get_nword() + 3) << " "
			<< *(datapoint + tko.get_nword() + 4) << " "
			<< *(datapoint + tko.get_nword() + 5) << " "
			<< std::endl;
		*/

		tko2.set_data(data + tko1.get_nword() + 5);
		ev_num2 = tko2.get_event_number(2);
		tko2_tag_event = ev_tag2 = tko2.get_event_tag(2);
		tko2_tag_spill = tko2.get_spill_tag(2);

		/*
		unsigned int * datapoint = data + tko1.get_nword() + 5;
		std::cout << "#D TKO2 " << std::hex
			<< *(datapoint + 0) << " "
			<< *(datapoint + 1) << " "
			<< *(datapoint + 2) << " "
			<< *(datapoint + 3) << " "
			<< *(datapoint + 4) << " "
			<< *(datapoint + 5) << " "
			<< std::endl;
		*/

		tko3.set_data(data + tko1.get_nword() + tko2.get_nword() + 5 + 4);
		//ev_num3 = tko3.get_event_number(2);
		//tko3_tag_event = ev_tag3 = tko3.get_event_tag(2);
		//tko3_tag_spill = tko3.get_spill_tag(2);

		/*
		datapoint = data + tko1.get_nword() + tko2.get_nword() + 5 + 4;
		std::cout << "#D TKO3 " << std::hex
			<< *(datapoint + 0) << " "
			<< *(datapoint + 1) << " "
			<< *(datapoint + 2) << " "
			<< *(datapoint + 3) << " "
			<< *(datapoint + 4) << " "
			<< *(datapoint + 5) << " "
			<< std::endl;
		*/

		if ((ev_num1 != ev_num2)
			|| ((ev_tag1 & 0xf) != (ev_tag2 & 0xf))) {
			std::cerr << "#E TKO event number miss match" << std::endl;
			std::cerr << "#E tko1 evnum:" << std::dec << ev_num1 << " ev_tag:" << ev_tag1
			<< " tko serial:"<< tko1.get_serial()
			<< " tko nword:" << tko1.get_nword()
			<< std::endl;
			std::cerr << "#E tko2 evnum:" << std::dec << ev_num2 << " ev_tag:" << ev_tag2
			<< " tko serial:" << tko2.get_serial()
			<< " tko nword:" << tko2.get_nword()
			<< std::endl;
		}



#if 0
		std::cout << "#D TKO 1" << std::endl;
		for (int i = 0 ; i < 32 ; i++) {
			std::cout << " " << tko1.get_data_val(20, i);
		}
		std::cout << std::endl;
		std::cout << "#D TKO 2" << std::endl;
		for (int i = 0 ; i < 32 ; i++) {
			std::cout << " " << tko2.get_data_val(20, i);
		}
		std::cout << std::endl;
#endif

	}

	return 0;
}
