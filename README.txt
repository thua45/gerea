Run Gerea on local computer

Step 1. installation from source

Go to the GEREA source folder.

Compile the source files using

gfortran -c 643.F

g++ -c main.cpp

g++ -c gerea.cpp

g++ -o gerea main.o gerea.o 643.o -lgfortran

An executable file named gerea will be generated

Step 2. Usage

To launch the program, open a terminal window, go to the folder where gerea executable file is located, and type:

./gerea in order to get the help information of the program

gerea -n network -e expression -s session

*network file

The file lists the links for the microRNA-(trascription factor)-(target gene) network or links dumped from the gredb database. The format should be tab delimited, The fist line should be "#db_type=2" where each following line is described by its "microRNA", "transcription factor" and "target gene" or "#db_type=1" where each following line is described by its "regulator", "target".

*expression file

The file lists the genes of interest and their expression in the transcriptome. The format should be tab delimited, The fist line should be "#data_type=2" where each following line is described by its "gene", "fold change", and "q value".

*Output file

An output file named session.ger.txt and session.links.txt. The first line consists of description, and the following lines de-scribe the microRNA name, values of a to i as well as the corrected P values (P0-P2). The network output file contains the links for each network that has been analyzed.

Step 3. Example run

Go to the GEREA folder. open a terminal window and run the program using

gerea -n Mouse_GEREDB.db.txt -e miR-155_CD4+.txt -f 1.8 -q 0.05 -s test_run -d out

Two out put files named "test_run.ger.txt" and "test_run.links.txt" will be generated int the out folder

An example run for the gredb-microRNA

gerea -n mouse-ge-mir-lcf.db.txt -e miR-155_CD4+.txt -f 1.8 -q 0.05 -s test_run2 -d out