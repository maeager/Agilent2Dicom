function nii = fdf2nii(dataset)
% nii = fdf2nii(dataset)

    [fse hdr] = readfdf([dataset '.img']);

    nii = make_nii(fse, hdr.voxelsize, [], 4, dataset);
    nii.hdr.hist.sform_code = 1;
    nii.hdr.hist.srow_x(1:3) = hdr.orientation(1,:);
    nii.hdr.hist.srow_y(1:3) = hdr.orientation(2,:);
    nii.hdr.hist.srow_z(1:3) = hdr.orientation(3,:);