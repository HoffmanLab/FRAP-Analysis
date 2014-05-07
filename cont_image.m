function cont_image(image)

figure;
colormap('gray')
imagesc(image)
hold on;
freezeColors
colormap('gray') % can change to be colored if desire
contour(image)
colorbar;

end