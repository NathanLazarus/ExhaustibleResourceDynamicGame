#!/bin/bash
donefile="Done.txt"

done=$(cat "$donefile")

if [[ $done -eq 1 ]]
then
  exit 0
else
  exit 222
fi