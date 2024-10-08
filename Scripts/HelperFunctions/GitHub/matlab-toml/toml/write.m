% WRITE serialize MATLAB data as TOML and write to file
%
%   WRITE('file.toml', struct('key', 5)) writes the text `key = 5` to
%   the file `file.toml`.
%
%   See also TOML.ENCODE, TOML.READ

function write(filename, matl_strct)
  fid = fopen(filename, 'w');
  fprintf(fid, '%s', encode(matl_strct));
  fclose(fid);
end