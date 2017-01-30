% for details on parameter choice etc., see the Word document

% amplitudes
IE = [30 20 55 45 36 40 20 30 25 30 15 25 35 44 40]/70*1e-9;
II = [90 60 110 100 54 50 52 60 60 70 60 50 60 40 30]/70*1e-9;


g02 = 1.8e-9;
g1 = 1.2e-9;
C = 60e-12;
C/(g02+g1)

Vp = 0;
Vm = -70e-3;

for factor = 1:500
w = factor*2*pi/10;

t = (1:5000)*1e-4;
ft = sin(t*w);

% figure, plot(t,ft)

dt = 1e-4;
for g = 2
    g3 = mean(IE)+(g-2)*std(IE);
    g4 = mean(II)+(g-2)*std(II);
    g3 = g3*0.5; % fudge factor >> 0.5 would be straight
    g4 = g4*0.5; % fudge factor >> 0.5 would be straight
    for kk = 1:3000
        k = (kk-1)*1e-3/100;
        V = -40e-3;
        for tt = 1:20000
            V = V + dt/C*( g02*(Vm-V) + g1*(Vp-V) + g3*(Vp-V)*sin(w*tt/1e4) + g4*(Vm-V)*sin(w*(tt/1e4 - k))  );
            VV(tt) = V;
        end
        amplitude(kk) = (max(VV(8000:20000))-min(VV(8000:20000))); 
    end
%     figure(12), plot((1:300)/10,amplitude*1e3); hold on;
    AA(:,factor) = amplitude;
end
end
hold off;
figure, imagesc([1:50],[1 300]/10,AA(:,:)*1000)
set(gca,'FontSize',13)
% figure(1), imagesc(u); colorbar
figure(2),hold on;
for g = 4
    plot(kkk,std(u(:,:,g)')*1000);
end
