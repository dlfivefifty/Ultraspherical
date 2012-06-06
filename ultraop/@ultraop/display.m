function display(N)
% DISPLAY Pretty-print a ultraop.
% DISPLAY is called automatically when a statement that results in a ultraop
% output is not terminated with a semicolon.

loose = ~isequal(get(0,'FormatSpacing'),'compact');
if loose, disp(' '), end
disp([inputname(1) ' = ultraop']);
if loose, disp(' '), end

fprintf('   Linear operator acting on functions defined on [-1,1]\n');
fprintf('   with n = 6 realisation:\n\n'); 
disp(realisation(N,6))

end