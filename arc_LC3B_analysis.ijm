// @File(label = "Input directory", style = "directory") inputdir
// @File(label = "Output directory", style = "directory") outputdir
// @String(label = "File suffix", value = ".tif") suffix

// Note: DO NOT DELETE OR MOVE THE FIRST 3 LINES -- they supply essential parameters

// arc_LC3B_analysis.ijm
// ImageJ macro to identify DAPI nuclei in somewhat faint/noisy images at high magnification,
// and measure nuclear and cytoplasmic intensity in two other channels
// Designed for analysis of neuronal markers
// Theresa Swayne, Ph.D, Columbia University, tcs6@cumc.columbia.edu

// Input: folder of single 3- channel images (usually max projections) with nuclei as the 1st channel
// Output: split channels, ROI set containing nuclei, measurements of nuclear and cytoplasmic intensity in the other channels
//
// Usage: Organize images to be analyzed into a folder. Run the macro.

// Limitations: only multi-channel images can be in the input directory. do not place previously split images there.
 
// ADJUSTABLE PARAMETERS -------------------------

// Neighborhood values -- larger nuclei require larger values
BLOCKSIZE = 127; // used in contrast enhancement
RADIUS = 40; // used in local thresholding

// Allowable nuclei sizes in microns^2
CELLMIN = 30; // minimum area
CELLMAX = 500; // maximum area

// Radius beyond the nucleus, in microns, that is used to measure cytoplasm
// Larger values may impinge on neighboring cells
// Smaller values may bring in more noise because of fewer pixels measured
CYTOPLASM_THICKNESS = 1

// --------- SETUP 

run("Set Measurements...", "area mean centroid integrated display decimal=2");
run("Clear Results");

// save data as csv, omit headers, preserve row number
run("Input/Output...", "file=.csv copy_row"); 

// add headers to results file
// 0 filename, 1 x centroid, 2 y centroid,
// 3-6 C2 nuclear, 7-10 C2 cyto,
// 11-14 C3 nuclear, 15-18 C3 cyto
headers1 = "Filename,X,Y,";
headers2 = "C2NucArea,C2NucMean,C2NucIntDen,C2NucRawIntDen,C2CytoArea,C2CytoMean,C2CytoIntDen,C2CytoRawIntDen,";
headers3 = "C3NucArea,C3NucMean,C3NucIntDen,C3NucRawIntDen,C3CytoArea,C3CytoMean,C3CytoIntDen,C3CytoRawIntDen";
headers = headers1 + headers2 + headers3;
File.append(headers,outputdir  + File.separator+ "Results.csv");

setBatchMode(true);
n = 0;

splitChannelsFolder(inputdir); // split each image into channels
processFolder(inputdir); // actually do the analysis
print("number of results at end of program:",nResults);
run("Clear Results");
print("number of results after final clearing:",nResults);

setBatchMode(false);

// ------- functions for processing folders

function splitChannelsFolder(inputdir) 
	{
   list = getFileList(inputdir);
   for (i=0; i<list.length; i++) 
   		{
   		isSingleChannel = ((startsWith(list[i], "C1")) || (startsWith(list[i], "C2")) || (startsWith(list[i], "C3")));
//   		print(list[i],isSingleChannel);	
        if(File.isDirectory(inputdir + File.separator + list[i])) {
			splitChannelsFolder("" + inputdir +File.separator+ list[i]);}
        else if (!isSingleChannel && (endsWith(list[i], suffix))) {     // avoids error if there are C1, C2, C3 images in the folder
			splitChannelsImage(inputdir, list[i]);}
    	}
	}

function processFolder(inputdir) 
	{
	list = getFileList(inputdir);
	for (i=0; i<list.length; i++) 
   		{
        if(File.isDirectory(inputdir + File.separator + list[i])){
			processFolder("" + inputdir +File.separator+ list[i]);}
        else if (endsWith(list[i], suffix))
        	{
        	if (startsWith(list[i], "C1")){
           		segmentNucleiImage(inputdir, list[i]);} // nuclei segmentation 
           	else if (startsWith(list[i], "C2")){
           		processC2C3Image(inputdir, list[i]);} // nuclei/cytoplasm intensity analysis
        	} // nothing happens with C3 images or original images
    	}
	}

// ------- functions for processing individual files

function splitChannelsImage(inputdir, name) 
	{
	open(inputdir+File.separator+name);
	print("splitting",n++, name);
	run("Split Channels");
	while (nImages > 0)  // works on any number of channels
		{
		saveAs ("tiff", inputdir+File.separator+getTitle);	// save every picture in the *input* folder
		close();
		}
	}

function segmentNucleiImage(inputdir, name) 
	{
	// assumes nuclei are in channel 1 of the previously split image, and the filename begins with "C1"
	open(inputdir+File.separator+name);
	print("opening C1 image",name);
	dotIndex = indexOf(name, ".");
	basename = substring(name, 3, dotIndex); // taking off the channel number
	procName = basename + "_processed.tif";
	resultName = basename + "_results.csv";
	nucRoiName = basename + "_Nuclei" + ".zip";
	id = getImageID();
	
	// process a copy of the image
	selectImage(id);
	// square brackets allow handing of filenames containing spaces
	run("Duplicate...", "title=" + "[" +procName+ "] duplicate"); 
	selectWindow(procName);
	
	// PRE-PROCESSING -----------------------------------------------------------
	run("Enhance Local Contrast (CLAHE)", "blocksize=" + BLOCKSIZE + " histogram=256 maximum=3 mask=*None*"); // accentuates faint nuclei
	run("Gaussian Blur...", "sigma=8"); // merges speckles to make nucleus more cohesive
	
	// SEGMENTATION AND MASK PROCESSING -------------------------------------------
	selectWindow(procName);
	run("Auto Local Threshold", "method=Phansalkar radius=" + RADIUS + " parameter_1=0 parameter_2=0 white");
	run("Convert to Mask");
	
	selectWindow(procName);
	run("Watershed"); // separate touching nuclei
	
	// analyze particles to get initial ROIs
	roiManager("reset");
	run("Analyze Particles...", "size=" + CELLMIN + "-" + CELLMAX + " exclude clear add");
	
	print("about to save the nuclear ROIs");// BUG: somewhere after this point the selection list is empty if C images pre-exist

	roiManager("Save", outputdir + File.separator + nucRoiName); 

	// clean up
	selectWindow(procName);
	close(); // BUG: it saves the processed image during a repeat run
	selectWindow(name);
	close();
	roiManager("reset");
	run("Clear Results");
	print("number of results at end of C1 function",nResults);
	}

function processC2C3Image(inputdir, name)
	{
	// converts nuclear ROIs to bands, 
	// measures nuclear and cytoplasmic (band) intensity in C2 and C3,
	// and saves results in a CSV file

	print("number of results at beginning of C2C3 function:",nResults);

	open(inputdir+File.separator+name); // C2 image
	print("opening C2 image",name);
	dotIndex = indexOf(name, ".");
	basename = substring(name, 3, dotIndex); // taking off the channel number
	resultName = basename + "_Results.csv";
	nucRoiName = basename + "_Nuclei.zip";
	cytRoiName = basename + "_Cyto.zip";

	C3Name = "C3-"+basename+".tif";
	open(inputdir+File.separator+C3Name); // corresponding C3 image
	print("opening C3 image",name);

	roiManager("Open", outputdir + File.separator + nucRoiName); // ROI set containing nuclei

	// measure C2 nuclear intensity
	selectWindow(name);
	roiManager("multi-measure measure_all append"); // measure individual nuclei and appends results
	// TODO: Is there a bug in multi-measure preventing clearing of results when append is selected?
	print("Measuring nuclei for",name);
	print("number of results after C2 nuclei:",nResults);

	// measure C3 nuclear intensity
	selectWindow(C3Name);
	roiManager("multi-measure measure_all append"); // measure individual nuclei and appends results
	print("Measuring nuclei for",C3Name);
	print("number of results after C3 nuclei:",nResults);
	
	// create and save cytoplasm band ROIs
	numROIs = roiManager("count");
	roiManager("Deselect");
	run("Select None");
	roiManager("Show None");
	for (index = 0; index < numROIs; index++) // loop through ROIs
		{
		roiManager("Select", index);
		run("Make Band...", "band="+CYTOPLASM_THICKNESS);
		roiManager("Update");
		}

	roiManager("Deselect");
	run("Select None");
	roiManager("Save", outputdir + File.separator + cytRoiName); 

	// measure C2 cytoplasmic intensity
	selectWindow(name);
	roiManager("multi-measure measure_all append");
	print("Measuring cytoplasm for",name);
	print("number of results after C2 cyto:",nResults);

	// measure C3 cytoplasmic intensity
	selectWindow(C3Name);
	roiManager("multi-measure measure_all append");
	print("Measuring cytoplasm for",C3Name);
	print("number of results after C3 cyto:",nResults);

	// TODO: write results

	// loop through numROIs (1 line per cell)
	
		// gather data
			// results table (n = roiManager(count)):
			// column 0 row, 1 label, 2 area, 3 mean, 4-5 x/y, 6-7 intden/rawintden
			// row 1 to n = C2 nuclei
			// n+1 to 2*n = C3 nuclei
			// (2*n)+1 to 3*n = C2 cyto
			// (3*n)+1 to 4*n = C3 cyto

		
		// assemble string resultString
			// 0 filename, 1 x centroid, 2 y centroid,
			// 3-6 C2 nuclear, 7-10 C2 cyto,
			// 11-14 C3 nuclear, 15-18 C3 cyto

		//  write data: 
			// File.append(resultString,outputdir + File.separator + resultName);

	


		// String.copyResults;
		// newResults = String.paste;
		// newResults = substring(newResults,0,lengthOf(newResults)-1); // strip the final newline
		// newResults = replace(newResults, "\t",","); // replace tabs with commas for csv

		// or getResultString("Column", row), getResult("Column", row), getResultLabel(row)
		// since using the actual results table

		//	selectWindow("Summary"); 
		//	lines = split(getInfo(), "\n"); 
		//	headings = split(lines[0], "\t"); 
		//	C1Values = split(lines[1], "\t"); 
		//	C1Name = C1Values[0];
		//	C2Values = split(lines[2], "\t"); 
		//	C2Name = C2Values[0];
		//	OverlapValues = split(lines[3], "\t"); 
		
			// convert strings to integers
		//	C1Count = parseInt(C1Values[1]);
		//	C2Count = parseInt(C2Values[1]);
		//	OverlapCount = parseInt(OverlapValues[1]);
		
			// calculate the percent overlap and convert to string
		//	C1withC2 = OverlapCount/C1Count;
		//	C2withC1 = OverlapCount/C2Count;
		//	strC1withC2 = d2s(C1withC2, 3);
		//	strC2withC1 = d2s(C2withC1, 3);
		
			// collect the colocalization data
			// format: image name, n cells, n cells colocalized with other label, % colocalized with other label
		//	firstRow = C1Name+",1," + C1Count+ ","+ OverlapCount + "," + strC1withC2;
		//	secondRow = C2Name + ",2," + C2Count+ "," + OverlapCount + "," + strC2withC1;
		//	colocResults = firstRow + "\n" + secondRow;


	// clean up
	selectWindow(name);
	close();
	selectWindow(C3Name);
	close();
	roiManager("reset");
	print("number of results at end of C2C3 function:",nResults);
	}

