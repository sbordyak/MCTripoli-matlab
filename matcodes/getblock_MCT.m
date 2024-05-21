function d = getblock_MCT(iblock,d0)




%d0.det(d0.det==1)=11;  %Move PM to back

[dets,~,newdets] = unique(d0.det);
isos = unique(d0.iso);

inb = d0.block==iblock;

d.data = d0.int(inb);
d.iso_vec = d0.iso(inb);
d.det_vec = newdets(inb);
d.blflag = ~d0.isOP(inb);
d.sig_ind = d.det_vec;
d.time = d0.time(inb);
d.axflag = d0.det(inb)==10; 
d.cycle = d0.cycle(inb); 

d.Ndata = length(d.data);
d.Ncycle = length(unique(d.cycle(d.cycle>0)));



d.Nfar = length(dets(dets~=10)); % 10 indicates PM
d.Niso = length(isos(isos>0));
d.Nblock = 1;
d.Ndet = length(dets);
d.Isotopes = unique(d0.mass(d0.iso~=0));
d.Nsig = size(d.sig_ind,2);



[timeblock,ininds,time_ind] = unique(d0.time(inb));
cycletmp = d.cycle(ininds);

time_ind(d.cycle==0) = 0;
mincycletime = min(time_ind(d.cycle~=0));
time_ind(d.cycle~=0) = time_ind(d.cycle~=0) + 1 - mincycletime;

d.time_ind = time_ind;




d.ReportInterval = mode(diff(timeblock));

%%

    for n = 1:d.Ncycle
        minCycleTime(n,1) = min(timeblock(cycletmp==n));
        maxCycleTime(n,1) = max(timeblock(cycletmp==n));
        
        iminCT(n,1) = find(timeblock==minCycleTime(n));
        imaxCT(n,1) = find(timeblock==maxCycleTime(n));
        
    end

    Tknots0 = [minCycleTime; maxCycleTime(end)];
    iTknots0 =  [iminCT; imaxCT(end)];






    %Block_Time{iblock} = Time(iTknots0(iblock,1):iTknots0(iblock,end));
    timeblock_noBL = timeblock(iTknots0(1):iTknots0(end));
    %InterpMat = interp1(Tknots0(iblock,:),eye(length(Tknots0(iblock,:))),Block_Time{iblock},'spline');
    InterpMat = interp1(Tknots0,eye(length(Tknots0)),timeblock_noBL,'linear');

    Nknots = length(Tknots0);
    Ntb = length(timeblock_noBL);

    %sb726 Added this to include cycle number for InterpMat/time index
    CycleMat = zeros(Ntb,1);
    for n = 1:d.Ncycle
        CycleMat((iminCT(n):imaxCT(n))-iminCT(1)+1) = n;
    end

    
%%



d.Nknots = Nknots;
d.Ntb = Ntb;

d.Include = true(d.Ndata,1);
d.IncludeMat{1}=true(d.Ntb,1);

d.InterpMat{1} = InterpMat;
d.CycleMat{1} = CycleMat;

d.block = true(size(d.data));


for ii = 1:d.Niso
    d.iso_ind(:,ii)=d.iso_vec==ii;
end
for ii = 1:d.Ndet
    d.det_ind(:,ii)=d.det_vec==ii;
end


d.Nt = length(unique(d.time));


end


