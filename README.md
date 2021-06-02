### MPCS-56420 Final Project

#### Edris Qarghah

#### 06/01/2021

##### How to Run

This project is an web visualization of the Needleman-Wunsch and Smith-Waterman sequence alignment algorithms. To run locally, first make sure you're in an environment running python3, then `pip install -r requirements.txt`. You can then simply run `main.py`.

You can also find the project online for the time being, hosted [here](https://seismic-hope-272623.uc.r.appspot.com). It may not be the most recent version of the app (i.e., it may have bugs not present if you run the provided code locally).

##### Functionality

Input any two sequences of arbtrary length. Once submitted, an empty table will populate and you can specify the algorithm, scoring matrix and gap cost. After doing so, you can step through the process of filling in the table, either a cell or a row at a time. If you select to Fill the entire table, you can then begin tracing the highest scoring alignment.

For Needleman-Wunsch (global alignment), this will start in the bottom right corner in the table and will proceed to the top right (arrows are not displayed on the for the first row (all come from the value to their left) and first column (all come from the cell above them). For Smith-Waterman, it will start at the maximum value across the entire table and trace back until it hits a zero.

##### Known funkiness

I appear to unintentionally be saving some details between runs. This means that if you change the sequence alignment, for example, your trace-back arrow may be wrong for a single frame. The indexing was working correctly before the class demo, so I'm guessing I messed something up in my last set of edits. Not sure if the version with the indexing issue is the one on the repo or not.

##### Post Deadline Update

I decided to go back and figure out what the indexing issue was, now that I have had some sleep. It was two-fold: I had not "fixed" the issue on all of my templates AND the orientation of seq1 and seq2 were reversed in creating the DP table (relative to what would have been expected from the input). Either way, it is fixed now.