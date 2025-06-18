// Enable batch processing mode to suppress GUI updates
setBatchMode(true);

// Define input and output directories
output = getDirectory("/Users/bowmanch/Desktop/CHP Auto/");
input  = getDirectory("/Users/bowmanch/Desktop/CHP Auto/Tiff test/");

// Ensure paths end with "/"
if (!endsWith(input, "/"))  input  += "/";
if (!endsWith(output, "/")) output += "/";

// Prompt user for file suffix
Dialog.create("File Type");
Dialog.addString("File Suffix: ", ".tiff");
suffix = Dialog.getString();

// Start processing
processFolder(input);

// === Function Definitions === //

function processFolder(input) {
    list = getFileList(input);
    for (i = 0; i < list.length; i++) {
        // Recursively process subfolders
        if (File.isDirectory(input + list[i])) {
            processFolder(input + list[i] + "/");
        }

        // Process files matching suffix
        if (endsWith(list[i], suffix)) {
            processFile(input, output, list[i]);
        }
    }
}

function processFile(input, output, file) {
    open(input + file);
    print("Processing: " + input + file);

    // Convert to RGB stack and apply colormap
    run("RGB Stack");
    run("mpl-viridis");

    // Measure basic properties
    run("Measure");
    run("Set Measurements...", "area mean standard min percent integrated limit display redirect=None decimal=3");

    // Apply threshold
    setAutoThreshold("Mean dark");
    setThreshold(50, 255, "raw");
    getThreshold(lower, upper);
    print("Min: " + lower + " - Max: " + upper);

    // Measure thresholded area
    run("Measure");

    // Save output as JPEG
    title = getTitle();
    saveAs("jpeg", output + title + ".jpg");
    print("Saved to: " + output + title + ".jpg");
}
