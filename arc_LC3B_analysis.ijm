// @File(label = "Input directory", style = "directory") dir1
// @File(label = "Output directory", style = "directory") dir2
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
 
// ADJUSTABLE PARAMETERS -------------------------

// Neighborhood values -- larger nuclei require larger values
BLOCKSIZE = 127; // used in contrast enhancement
RADIUS = 40; // used in local thresholding

// Allowable nuclei sizes in microns^2
CELLMIN = 30; // minimum area
CELLMAX = 500; // maximum area

// Radius beyond the nucleus, in microns, that is used to measure cytoplasm
// Larger values may impinge on neighboring cells
// Smaller values have more noise
CYTOPLASM_THICKNESS = 5

// SETUP -----------------------------------------------------------------------

run("Set Measurements...", "area mean min centroid integrated display decimal=2");

// save data as csv, preserve headers, preserve row number for copy/paste
run("Input/Output...", "file=.csv copy_row save_column save_row"); 

// add headers to results file
headers = "Filename,Channel,Nuc Area,Nuc Mean,Nuc IntDen,Nuc RawIntDen,Cyto Area,Cyto Mean,Cyto IntDen,Cyto RawIntDen";
File.append(headers,dir2  + File.separator+ "Results.csv");

setBatchMode(true);
n = 0;

splitChannelsFolder(dir1); // split each image into channels
processFolder(dir1); // actually do the analysis

setBatchMode(false);

function splitChannelsFolder(dir1) 
	{
   list = getFileList(dir1);
   for (i=0; i<list.length; i++) 
   		{
        if(File.isDirectory(dir1 + File.separator + list[i])) {
			processFolder("" + dir1 +File.separator+ list[i]);}
        else if (endsWith(list[i], suffix)) {
           		splitChannelsImage(dir1, list[i]);}
    	}
	}


function processFolder(dir1) 
	{
   list = getFileList(dir1);
   for (i=0; i<list.length; i++) 
   		{
        if(File.isDirectory(dir1 + File.separator + list[i])){
			processFolder("" + dir1 +File.separator+ list[i]);}
        else if (endsWith(list[i], suffix))
        	{
        	if (startsWith(list[i], "C1")){
           		segmentNucleiImage(dir1, list[i]);} // gets the nuclei -- assumes filename begins with "C1"
           	else if (startsWith(list[i], "C2")){
           		processC2Image(dir1, list[i]);} // TODO: nuclei/cytoplasm analysis
        	}
    	}
	}

function splitChannelsImage(dir1, name) 
	{
	open(dir1+File.separator+name);
	print("splitting",n++, name);
	run("Split Channels");
	while (nImages > 0)  // works on any number of channels
		{
		saveAs ("tiff", dir1+File.separator+getTitle);	// save every picture in the *input* folder
		close();
		}
	}

function segmentNucleiImage(dir1, name) 
	{
	// assumes nuclei are in channel 1 of the previously split image, and the filename begins with "C1"
	open(dir1+File.separator+name);
	dotIndex = indexOf(name, ".");
	basename = substring(name, 3, dotIndex); // taking off the channel number
	procName = "processed_" + basename + ".tif";
	resultName = "results_" + basename + ".csv";
	roiName = "Nuclei_" + basename + ".zip";
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
	//run("Options...", "iterations=" + OPENITER + " count=" + OPENCOUNT + " black"); // smooth borders
	//run("Open");
	run("Watershed"); // separate touching nuclei
	
	// analyze particles to get initial ROIs
	
	roiManager("reset");
	run("Analyze Particles...", "size=" + CELLMIN + "-" + CELLMAX + " exclude clear add");
	roiManager("Save", dir2 + File.separator + roiName); // saved in the output folder
	}

function processC2Image(dir1, name)
	{
	open(dir1+File.separator+name);
	dotIndex = indexOf(name, ".");
	basename = substring(name, 3, dotIndex); // taking off the channel number
	procName = "processed_" + basename + ".tif";
	resultName = "results_" + basename + ".csv";
	roiName = "Nuclei_" + basename + ".zip";
	id = getImageID();

	roiManager("Open", dir1 + File.separator + roiName); // open the ROI set containing nuclei

	// measure nuclear intensity
	roiManager("multi-measure measure_all append"); // measures individual nuclei and appends results

	// TODO: measure cytoplasmic intensity
		// generate band ROIs

		// update and save cytoplasmic rois

		// measure intensity

	// write results by cannibalizing the below

		//	String.copyResults;
		//	newResults = String.paste;
		//	newResults = substring(newResults,0,lengthOf(newResults)-1); // strip the final newline
		//	newResults = replace(newResults, "\t",","); // replace tabs with commas for csv

		
				// read data from the Summary window
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

		//	File.append(newResults,dir2 + File.separator + resultName);
			
	}

