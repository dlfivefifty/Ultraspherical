function rows = BCrows(Wa,Wb,n)
% returns entries in the matrix required to impose the robin boundary rows.


if size(Wa,2) > 2
    error('BOUNDARY:ROWS','Currently do not support high order boundary conditions.');
end
if size(Wa,1) > 2
    error('BOUNDARY:ROWS','Currently only support boundary conditions at endpoints.');
end

rows = zeros(size(Wa,1),n);

% template
leftdirichlet = (-1).^(0:n-1);
rightdirichlet = ones(1,n);
leftneumann = (-1).^(1:n).*(0:n-1).^2;
rightneumann = (0:n-1).^2;

% dirichlet conditions.
rows = rows + [leftdirichlet;leftdirichlet].*kron(Wa(:,1),ones(1,n));
rows = rows + [rightdirichlet;rightdirichlet].*kron(Wb(:,1),ones(1,n));

% neumann condtions.
if size(Wa,2) > 1
    rows = rows + [leftneumann;leftneumann].*kron(Wa(:,2),ones(1,n));
    rows = rows + [rightneumann;rightneumann].*kron(Wb(:,2),ones(1,n));
end

end