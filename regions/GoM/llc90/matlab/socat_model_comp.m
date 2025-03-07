addpath('C:\Users\jessicaz\Documents\MATLAB\lib')

load gom_grid_data.mat


nx=20;
ny=15;

dataDir = 'C:\Users\jessicaz\Documents\research\regions\GoM\llc90\fugfCO2\'; 

mod_fco2=nan(nx,ny);

for rr=1:length(R)

    year=R(rr);

    for tmpMon=1:12

        fnpre=['fugfCO2_' num2str(nx) 'x' num2str(ny) '..'];
        fnsuf=['T120000'];
        
        if tmpMon<10
            Files = dir([dataDir fnpre num2str(year) '0' ...
                num2str(tmpMon) '*' fnsuf]);
        else 
             Files = dir([dataDir fnpre num2str(year) ...
                 num2str(tmpMon) '*' fnsuf]);
        end
  
    
        mdata=readbin([dataDir Files.name],[nx ny]);
    
        mod_fco2=cat(3,mod_fco2,mdata);

    end
end


figure()
scatter(gom_fco2(:),mod_fco2(:))






