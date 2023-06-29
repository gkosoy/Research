function action(input, output, filename) {
        open(input + filename);
        makeRectangle(70, 60, 630, 700);
        run("Crop");
        saveAs("tiff", output + filename);
        close();
}

path1 = "C:/Users/15183/Documents/summer_2023/tocrop/cal_ser_3_20.tiff";
path2 = "C:/Users/15183/Documents/summer_2023/tocrop/cropped/new.tiff";
input = File.getDirectory(path1);
output = File.getDirectory(path2);

setBatchMode(true); 
list = getFileList(input);
for (i = 0; i < list.length; i++){
        action(input, output, list[i]);
}
setBatchMode(false);





