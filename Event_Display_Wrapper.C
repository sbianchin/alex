
ifstream fdat_read("/Users/alng/Documents/Datascience/TREK/offline/alex/Keito_run_3994_1_test.txt", ios::in);


/*
*/

int Event_Display_Wrapper () {
  while (fdat_read.good())
  {
    int Read_Run_Number;
    fdat_read >> Read_Run_Number;
    Event_Display_5_1(3994,Read_Run_Number);
  }
  return 0;
}
