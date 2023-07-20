function Q = qcritcalc(hx,hy,hz,u)
 
Q = 0*squeeze(u(:,:,:,1,:));
 
for tt = 1:size(u,5)
     
    [dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz] ...
        = derivatives3d(hx, hy, hz, u(:,:,:,:,tt));
     
    S{1,1} = 0.5 * (dudx + dudx);
    S{1,2} = 0.5 * (dudy + dvdx);
    S{1,3} = 0.5 * (dudz + dwdx);
     
    S{2,1} = 0.5 * (dvdx + dudy);
    S{2,2} = 0.5 * (dvdy + dvdy);
    S{2,3} = 0.5 * (dvdz + dwdy);
     
    S{3,1} = 0.5 * (dwdx + dudz);
    S{3,2} = 0.5 * (dwdy + dvdz);
    S{3,3} = 0.5 * (dwdz + dwdz);
     
    O{1,1} = 0.5 * (dudx - dudx);
    O{1,2} = 0.5 * (dudy - dvdx);
    O{1,3} = 0.5 * (dudz - dwdx);
     
    O{2,1} = 0.5 * (dvdx - dudy);
    O{2,2} = 0.5 * (dvdy - dvdy);
    O{2,3} = 0.5 * (dvdz - dwdy);
     
    O{3,1} = 0.5 * (dwdx - dudz);
    O{3,2} = 0.5 * (dwdy - dvdz);
    O{3,3} = 0.5 * (dwdz - dwdz);
     
    Snorm = 0;
    Onorm = 0;
     
    for iii = 1:3
        for jjj = 1:3
             
            Snorm = Snorm + S{iii,jjj}.^2;
            Onorm = Onorm + O{iii,jjj}.^2;
             
        end
    end
     
    %Snorm = sqrt(Snorm);
    %Onorm = sqrt(Onorm);
     
    Q(:,:,:,tt) = 0.5*(Onorm - Snorm);
     
end