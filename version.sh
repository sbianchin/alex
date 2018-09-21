#!/bin/bash

rm -f Event_Display_C.d
rm -f Event_Display_C.so
rm -f Event_Display.C

ln -s Event_Display_$1.C Event_Display.C