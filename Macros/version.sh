#!/bin/bash

rm -f Event_DisplayC.d
rm -f Event_Display.so
rm -f Event_Display.C

ln -s Event_Display_$1.C Event_Display.C