# AutoMitoNetwork
See AutoMitoNetwork manual for full document on installation and instructions

Installation instructions
•	Open MyAppInstaller_web and follow prompts to download the application. 
•	Once completed, a folder labelled AutoMitoNetwork will appear in your selected save location. 
•	Open application folder.
•	Select AutoMitoNetwork to open application and start analysis.
 

Application analysis instructions

1.	Press load images to upload image(s) for analysis.
The number of images selected and the folder path they are in will appear. 
 
2.	Input the following optional parameters: 
Scale (1=default). If scale is known can input value here, otherwise units are in pixels. 
Allowable mitochondria size (pixels) (1=default). This removes detected mitochondrial networks with an area less than the inputted value. 

3.	Select shape and gray-level intensity features.
A checked box in the parent node ‘Shape features’ or ‘Gray-level intensity features’ means all features are selected. 
A coloured box in these parent nodes mean only certain features were selected.

4.	Check texture features if texture features are to be analysed. If checked proceed to the following step otherwise skip to step 8.
 
5.	Input offset for the gray-level co-occurrence matrices as a matrix. The offset is the distance between the pixel of interest and its neighbour. Each row in the matrix is a two-element vector, [row_offset, col_offset], that specifies the relationship, or offset, of a pair of pixels. row_offset is the number of rows between the pixel-of-interest and its neighbour. col_offset is the number of columns between the pixel-of-interest and its neighbour. See https://au.mathworks.com/help/images/ref/graycomatrix.html for more information. 
Default = [0 1; -1 1; -1 0; -1 -1], which is four offsets: one pixel to the right, one pixel up, one pixel 45° to the right and one pixel 45° to the left. The figure below shows how offsets can be defined. 
 
6.	Input number of gray levels for the gray-level co-occurrence matrices.
 
7.	Select texture features to be generated. 
 
8.	Click calculate to perform analysis and generation of features.
- A .xls file is generated and saved in the folder path containing the features for each detected mitochondria network per image. The name of the file is called filename_ _features_per_mitonetwork.xls
- A .png file is generated and saved in the folder path containing an image of the detected mitochondrial networks.  The name of the file is called filename_mitomap.png

9.	The figure panel in the app will show the final image that was analysed, showing the detected mitochondrial networks and that the program has finished calculations. 

