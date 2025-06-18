// Prompt user to select the directory containing TIFF images
path = getDirectory("/Users/bowmanch/Desktop/Exported TIFFS/");
filelist = getFileList(path);

// Clear any existing results and set batch mode
run("Clear Results");
setBatchMode(true);

// Loop through all files in the selected directory
for (i = 0; i < filelist.length; i++) {
    print("Analyzing: " + filelist[i]);

    // Process only TIFF files
    if (endsWith(filelist[i], ".tif") || endsWith(filelist[i], ".tiff")) {
        print("Filename is: " + path + filelist[i]);

        // Open the image
        open(path + filelist[i]);

        // Run color deconvolution using H DAB vector
        print("Color filename is: " + path + filelist[i] + "-(Colour_2)");
        run("Colour Deconvolution", "vectors=[H DAB]");

        // Select the deconvolved DAB image
        selectWindow(filelist[i] + "-(Colour_2)");

        // Set and apply threshold
        setAutoThreshold("Default");
        setThreshold(0, 192);
        setOption("BlackBackground", true);
        run("Convert to Mask");

        // Measure and store results
        run("Measure");
        print("Running it for filename: " + filelist[i]);
        setResult("Label", i, filelist[i]);

        // Close all open images
        run("Close All");
    }
}

// Print completion message
print("Ran all");

// Set output file name (results are automatically saved when using File > Save As in the GUI)
resultName = "results.xls";
