# AngioMT
A MATLAB based 2D image-to-physics tool to predict oxygen transport in vascularized microphysiological systems.

<h2>Dependencies</h2>
This program uses Im2mesh, an open source MATLAB script that creates finite element meshes from 2D grayscale images [1]. Please download this package and add to the working path for running AngioMT. <br><br>

<ol>
  <li>Jiexian Ma (2024). Im2mesh (2D image to triangular meshes) (https://www.mathworks.com/matlabcentral/fileexchange/71772-im2mesh-2d-image-to-triangular-meshes), MATLAB Central File Exchange.</li>
</ol>

<h2>Running the code</h2>
<ol>
  <li>Download the main program, <i>AngioMT_Main.m</i>, along with its dependencies (functions) <i>img_bin.m, element_quality.m, global_matrices.m and aux_var.m</i>.</li>
  <li>Open <i>AngioMT_Main.m</i> in MATLAB and add values of oxygen mass transport parameters (for ex. diffusivity, flux, reaction constant etc.) in lines 7 through 13 of <i>AngioMT_Main.m</i>.</li>
  <li>Download <i>Im2mesh</i> and add the path to its folder in line 16 of <i>AngioMT_Main.m</i>.
  <li>Add the path to directory containing images to be analyzed in line 19 of <i>AngioMT_Main.m</i>.<br><br><b>NOTE:</b> The current version has been developed to read '.tif' images. For other image types, please change the extension of the filename in line 20.<br><b>NOTE:</b> If users have pre-binarized images from other packages (REAVER etc.), remove the image binarization step in line 39 of <i>AngioMT_Main.m</i>.<br><br>
  <li>Run the program by clicking "Run" command in MATLAB.</li>
</ol>

<h2>Analyzing results</h2>
By default, AngioMT generates three figures:
<ol>
  <li> <b>Image processing:</b> Original image, binarized image and the segmented image with newe grayscale values assigned to disconnected domains (ex. disconnected vessels).</li>
  <li> <b>Meshing:</b> Meshed networks after Delaunay triangulation and the element quality map.</li>
  <li> <b>Mass transport results:</b> Computed oxygen distribution and oxygen flux plots.</li>
</ol>

AngioMT also creates a global table storing ancillary variables (Vascular oxygen, Oxygen delivery, Flux) for each image analyzed during the routine.
