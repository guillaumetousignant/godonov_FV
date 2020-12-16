function make_su2_mesh(varargin)
    %MAKE_SU2_MESH This program creates a square mesh.
    %   The input is made with parameter-value pairs, using the following syntax: 
    %   "make_su2_mesh('filename', 'meshes/mesh.su2', 'xmin', 0.5, 'PARAMETER', VALUE, [...]);"
    %   All parameters are optional, defaults are set otherwise.
    %
    %       Parameter   [default value]         Description 
    %
    %       filename    ['meshes/mesh.su2']     Output filename or path to mesh file to be written.
    %       xmin        [-0.5]                  Mesh span in the x direction start.
    %       xmax        [0.5]                   Mesh span in the x direction end.
    %       ymin        [-0.5]                  Mesh span in the y direction start.
    %       ymax        [0.5]                   Mesh span in the y direction end.
    %       xres        [500]                   Number of cells in the x direction
    %       yres        [500]                   Number of cells in the y direction
    %       boundary    [farfield]              Tag to use for the boundaries

        % Default values
        filename = 'meshes/mesh.su2';
        xmin = -0.5;
        xmax = 0.5;
        ymin = -0.5;
        ymax = 0.5;
        xres = 500;
        yres = 500;
        boundary_tag = 'farfield';
    
        % Input parsing
        if ~isempty(varargin)
            if rem(length(varargin), 2)
                error('make_circle_su2:unevenArgumentCount', 'Error, uneven argument count. Arguments should follow the "''key'', value" format. Exiting.');
            end
            for i = 1:2:length(varargin)
                key = varargin{i};
                value = varargin{i+1};
    
                switch lower(key)
                    case "filename"
                        filename = value;
                    case "xmin"
                        xmin = value;
                    case "xmax"
                        xmax = value;
                    case "ymin"
                        ymin = value;
                    case "ymax"
                        ymax = value;
                    case "xres"
                        xres = value;
                    case "yres"
                        yres = value;
                    case "boundary"
                        boundary_tag = value;
                    otherwise
                        warning('Warning, unknown parameter: ''%s'', ignoring.', key);
                end
            end
        end

        n_elements = xres * yres;
        n_points = (xres + 1) * (yres + 1);
        n_boundaries = 2 * xres + 2 * yres;

        x = linspace(xmin, xmax, xres + 1);
        y = linspace(ymax, ymin, yres + 1); % Backwards because I made my drawing backwards and I don't want to draw it again.

        points = zeros(n_points, 2);
        elements = zeros(n_elements, 4);
        boundaries = zeros(n_boundaries, 2);

        for i = 1:xres + 1
            for j = 1:yres + 1
                points(i + (j - 1) * (xres + 1), 1) = x(i);
                points(i + (j - 1) * (xres + 1), 2) = y(j);
            end
        end

        for i = 1:xres
            for j = 1:yres
                elements(i + (j - 1) * xres, 1) = i + (j - 1) * (xres + 1);
                elements(i + (j - 1) * xres, 4) = i + (j - 1) * (xres + 1) + 1;
                elements(i + (j - 1) * xres, 2) = i + j * (xres + 1);
                elements(i + (j - 1) * xres, 3) = i + j * (xres + 1) + 1;
            end
        end

        for i = 1:xres
            boundaries(i, 2) = i;
            boundaries(i, 1) = i + 1;

            offset = xres;
            boundaries(i + offset, 1) = yres * (xres + 1) + i;
            boundaries(i + offset, 2) = yres * (xres + 1) + i + 1;
        end

        for j = 1:yres
            offset = 2 * xres;
            boundaries(j + offset, 2) = j * (xres + 1);
            boundaries(j + offset, 1) = (j + 1) * (xres + 1);

            offset = 2 * xres + yres;
            boundaries(j + offset, 1) = (j - 1) * (xres + 1) + 1;
            boundaries(j + offset, 2) = j * (xres + 1) + 1;
        end
    
        write_su2(filename);
    
        function write_su2(filename)
            filepath = fileparts(filename);
            if ~isempty(filepath) && ~exist(filepath, 'dir')
                mkdir(filepath);
            end

            su2_file = fopen(filename, 'w');
        
            fprintf(su2_file, 'NDIME= 2\n\n');
            fprintf(su2_file, 'NPOIN= %g\n', n_points);
            for k = 1:n_points
                fprintf(su2_file, '%g %g\n', points(k, 1), points(k, 2));
            end
        
            fprintf(su2_file, '\nNELEM= %g\n', n_elements);
            for k = 1:n_elements
                fprintf(su2_file, '9 %g %g %g %g\n', elements(k, 1), elements(k, 2), elements(k, 3), elements(k, 4));
            end
    
            fprintf(su2_file, 'NMARK= 1\n');
            fprintf(su2_file, 'MARKER_TAG= %s\n', boundary_tag);
            fprintf(su2_file, 'MARKER_ELEMS= %g\n', n_boundaries);
            for k = 1:n_boundaries
                fprintf(su2_file, '3 %g %g\n', boundaries(k, 1), boundaries(k, 2));
            end
    
            fclose(su2_file);
        end
end
    