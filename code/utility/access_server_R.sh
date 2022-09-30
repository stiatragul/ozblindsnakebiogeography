#!/bin/bash
sshpass -f "/home/st_wsl/kak.txt" ssh putter@nick.rsb.anu.edu.au 'cd blindsnakebiogeography; R --vanilla'

# then run source'code/....' in the R console