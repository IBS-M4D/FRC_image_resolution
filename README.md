# FRC_image_resolution
Calculate SMLM image resolution with the Fourier Ring Correlation (FRC) algorithm by comparing spatial frequencies between two independent reconstructions. Works with ThunderSTORM results tables. 

How to use it:
1. Run the .exe file
2. Presse the Browse CSV button to search and load the ThunderSTORM results table 
The code splits the data in two halves to reconstruct two almost identical images, only one of them is displayed on the screen.
3. Draw a rectangular ROI on the reconstructed image
4. You can change the size of ROI and/or move it within the image to estimate the local resolution.
5. Pixel Size of the reconstructed image might bottleneck the estimation of the real resolution. You can decrease it form 10 to 5 nm (or less) to see if this improves the calculated resolution.

If you get a noisy curve or an error, it is oten related to the selected region or its size (not sufficient amount of data). Try to increase the ROI size, or select the whole image.
