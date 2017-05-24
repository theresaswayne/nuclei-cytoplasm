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
CYTOPLASM_THICKNESS = 1

// --------- SETUP 

run("Set Measurements...", "area mean min centroid integrated display decimal=2");
run("Clear Results");

// save data as csv, preserve headers, preserve row number for copy/paste
run("Input/Output...", "file=.csv copy_row save_column save_row"); 

// add headers to results file
headers = "Filename,Channel,Nuc Area,Nuc Mean,Nuc IntDen,Nuc RawIntDen,Cyto Area,Cyto Mean,Cyto IntDen,Cyto RawIntDen";
File.append(headers,outputdir  + File.separator+ "Results.csv");

setBatchMode(true);
n = 0;

splitChannelsFolder(inputdir); // split each image into channels
processFolder(inputdir); // actually do the analysis

setBatchMode(false);

// ------- functions for processing folders

function splitChannelsFolder(inputdir) 
	{
   list = getFileList(inputdir);
   for (i=0; i<list.length; i++) 
   		{
        if(File.isDirectory(inputdir + File.separator + list[i])) {
			processFolder("" + inputdir +File.separator+ list[i]);}
        else if (endsWith(list[i], suffix)) {
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
           		segmentNucleiImage(inputdir, list[i]);} // gets the nuclei -- assumes filename begins with "C1"
           	else if (startsWith(list[i], "C2")){
           		processC2Image(inputdir, list[i]);} // nuclei/cytoplasm analysis
        	}
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
	roiManager("Save", outputdir + File.separator + nucRoiName);

	// clean up
	selectWindow(procName);
	close();
	selectWindow(name);
	close();
	roiManager("reset");
	}

function processC2Image(inputdir, name)
	{
	// converts nuclear ROIs to bands and measures intensity

	open(inputdir+File.separator+name);
	dotIndex = indexOf(name, ".");
	basename = substring(name, 3, dotIndex); // taking off the channel number
	resultName = basename + "_Results.csv";
	nucRoiName = basename + "_Nuclei.zip";
	cytRoiName = basename + "_Cyto.zip";
	id = getImageID();

	roiManager("Open", outputdir + File.separator + nucRoiName); // open the ROI set containing nuclei

	// measure nuclear intensity
	roiManager("multi-measure measure_all append"); // measures individual nuclei and appends results

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

	// measure cytoplasmic intensity
	roiManager("multi-measure measure_all append"); // measures individual cytoplasms and appends results

	// save cytoplasm ROIs
	roiManager("Save", outputdir + File.separator + cytRoiName); // saved in the output folder

	// TODO: write results by cannibalizing the below

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

		//	File.append(newResults,outputdir + File.separator + resultName);

	// clean up
	selectWindow(name);
	close();
	roiManager("reset");
	}

