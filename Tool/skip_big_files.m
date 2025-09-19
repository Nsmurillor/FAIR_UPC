folder = 'D:\Github\FAIR_UPC';  % cambia por tu ruta

files = dir(fullfile(folder, '**', '*'));
files = files(~[files.isdir]);  % eliminar carpetas

% Tamaño en MB
fileSizesMB = [files.bytes] / 1e6;

% Filtrar archivos mayores a 100 MB
largeFiles = files(fileSizesMB > 60);

% Crear o abrir .gitignore
fid = fopen('big_files.txt', 'w');

for k = 1:length(largeFiles)
    % Path relativo y con /
    relPath = strrep(fullfile(largeFiles(k).folder, largeFiles(k).name), folder, '');
    relPath = strrep(relPath, '\', '/');   % reemplazar \ por /
    % Asegurarse que empieza sin /
    if startsWith(relPath, '/')
        relPath = relPath(2:end);
    end
    fprintf(fid, '%s\n', relPath);  % agregar / al inicio
end

fclose(fid)