array_num = 24;
d_min = 0.02;

cell0 = (1:(array_num-1)/3)*d_min;
cell3 = (1j)*cell0*exp(1j*pi)*exp(1j*pi/2);
cell2 = cell0*exp(-1j*pi/6)*exp(1j*pi)*exp(1j*pi/2);
cell1 = cell0*exp(-1j*pi*5/6)*exp(1j*pi)*exp(1j*pi/2);
ant_pos_x = [0 real([cell1 cell2 cell3])];
ant_pos_y = [0 imag([cell1 cell2 cell3])];
ant_pos = [ant_pos_x;ant_pos_y];
elem_pos_opt = ant_pos(1,:)+1j*ant_pos(2,:);
elem_pos_opt = elem_pos_opt.';
save('array_Y.mat','elem_pos_opt');
figure
scatter(real(elem_pos_opt),imag(elem_pos_opt));