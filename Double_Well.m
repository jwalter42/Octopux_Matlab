%*********************************************Double Well VNA ER Updated***************************************************
%test_flag=1, do ER only; flag=2, do VNA only; flag=3, do VNA then ER; flag=4, do ER then VNA

function [Data1, AvgPwr1, e1, e1_pat, e1_ns, e1_ds, e1_os, Data2, AvgPwr2, e2, e2_pat, e2_ns, e2_ds, e2_os]= Double_Well(FFL_or_Phip, phin_path1, phin_sweep1, phip_levels, power, avg, bw, fstart, fstop, npoints, sample_range, pattest, ntest, edges, bitwise, test_flag, tag, loop_times)   
    global imb iob icb itb idb ibs ffl
    global PATTERN pat1 pat1b pat2b pat3b pat4b pat100b pat10b pat400b pat1000b pat4000b pat10000b pat40000b pcl
    global imb_chan iob_chan itb_chan idb_chan icb_chan ibs_chan ffl_chan
    global v gnd fridge chassis
    %Sweep parameters
    if nargin<1 || isempty(FFL_or_Phip),  FFL_or_Phip=  "Phip"; end
    if nargin<2 || isempty(phin_path1),  phin_path1=  0:.01:0.8; end
    if nargin<3 || isempty(phin_sweep1), phin_sweep1= 0.8:-.01:-0.8; end
    if nargin<4 || isempty(phip_levels), phip_levels=  [0,0.4]; end
    %VNA parameters
    if nargin<5 || isempty(power), power= -45; end
    if nargin<6 || isempty(avg), avg= 2; end
    if nargin<7 || isempty(bw), bw=1e3; end
    if nargin<8 || isempty(fstart), fstart= 6.7e9; end
    if nargin<9 || isempty(fstop), fstop=6.7000001e9; end
    if nargin<10 || isempty(npoints), npoints= 31; end
    if nargin<11 || isempty(sample_range), sample_range=[fstart, fstop]; end
    %Comp parameters
    if nargin<12 || isempty(pattest), pattest=1; end
    if nargin<13 || isempty(ntest), ntest=1; end
    if nargin<14 || isempty(edges), edges=1; end
    if nargin<15 || isempty(bitwise), bitwise=true; end
    %Label/Flag
    if nargin<16 || isempty(test_flag), test_flag=2; end
    if nargin<17 || isempty(tag), tag=""; end
    if nargin<18 || isempty(loop_times), loop_times=1; end

    %Global equipment definitions/ error catches
    if ~isfinite(loop_times) || loop_times < 1 || loop_times ~= fix(loop_times)
        disp("loop_times not a positive integer, check variable definitions")
        return
    end
    if test_flag~=1 && (~exist('chassis', 'var') || isempty(chassis))
        chassis=input("Enter chassis number of Keysight VNA (usually 3 for ADRs): ");
    end
    if ~exist('fridge', 'var') || isempty(fridge)
        fridge=input("Enter fridge in use, for temp grabs ('adr1', 'adr2', 'He', 'dr', etc): ");
    end
    if any(sample_range < 1e6, 'all') || any([fstart, fstop] < 1e6, 'all')
        disp("Frequencies may not be defined correctly, less than 1 MHz");
    end

    %Comp initializing
    e1=[];
    e2=[];
    e1_pat=[];
    e2_pat=[];
    e1_ns=[]; 
    e1_ds=[];
    e1_os=[];
    e2_ns=[];
    e2_ds=[];
    e2_os=[];
    %splitcurrent = 0; %May be needed in future edits

    if test_flag~=2
        bits=sum(PATTERN(pattest).vdat,"all");
        if isnumeric(bitwise) && bitwise>0
            bits=round(bitwise);
        end
     end

    %VNA initializing
    f1={};
    Data1={};
    AvgPwr1={};
    f2={};
    Data2={};
    AvgPwr2={};   
   
    if test_flag ~= 1
        clear v
        v = PXI_VNA();
        v.setOutput(1);
        v.setSparam("S21");
        v.setPower(power)
        v.setAverage(avg)
        v.setBandwidth(bw)
        v.setStopFrequency(40e9)
        v.setStartFrequency(fstart)
        v.setStopFrequency(fstop)
        v.setNpts(npoints)
        v.getOutput
        v.getSparam
        v.getPower
        v.getAverage
        v.getBandwidth
        v.getStartFrequency
        v.getStopFrequency
        v.getNpts
        
        rang=linspace(fstart,fstop,npoints);  %Live results for first sample range only (if there are more than one)
        mini=sample_range(1,1);
        maxi=sample_range(1,2);
        idx1=find(rang>=mini,1,"first");
        idx2=find(rang<=maxi,1,"last");
        if isempty(idx1), idx1=1; end
        if isempty(idx2), idx2=length(rang); end
    end

    % icb is array of [DCFlux1, DCFlux2], also called icb(1) and icb(2), they are currents applied to induce flux in JJ1 and JJ2. 
    % Phi+= DCFlux1 - DCFlux2, Phi-= DCFlux1 + DCFlux2
    % DCFlux1 = (Phi- + Phi+)/2, DCFlux2 = (Phi- - Phi+)/2

    phin0=phin_path1(1);                    %Starting Phi minus for both sweeps
    phin_path2=2*phin_path1(1)-phin_path1;  %Define return initialization sweep based on first initialization sweep
    phin_sweep2=flip(phin_sweep1);          %Define return test sweep based on first test sweep
    
    if length(icb_chan)==1                  %With 1 icb_chan bias, assume 1-line Phi-. Also use FFl 
        Phin=1;
        FFL_or_Phip="FFL";
        Phin_string=string("1-line Phi-");
    elseif length(icb_chan)==2              %With 2 icb_chan bias, assume 2-line Phi- = DCFLux1 + DCFlux2. Can be used with Phi+= DCFLux1 - DCFlux2 or with ffl_chan
        Phin=0;
        Phin_string=string("2-line Phi-");
    end

    %Set test to use Phip or FFL; if using fast flux, redefine terms to use ffl0 and ffl_PI and set all phip variables to 0
    if istrue(FFL_or_Phip) || contains(FFL_or_Phip, "f", "IgnoreCase", true)
        FF=1;                            
        ffl0=phip_levels(1); 
        ffl_PI=phip_levels(2);
        phip0=0;
        phip_PI=0;
        PI_string=string(ffl_PI);
    elseif isfalse(FFL_or_Phip) || contains(FFL_or_Phip, "p", "IgnoreCase", true)
        FF=0;
        phip0=phip_levels(1);   %Phi plus corresponding to 0 state
        phip_PI=phip_levels(2); %Phi plus corresponding to PI state
        PI_string=string(phip_PI);
    else
        disp("Specify 'Phip' or 'FFL' as first variable")
        return
    end

    for loops =1:loop_times
        fprintf("*****************************Iteration %d Sweep1*******************\n",loops)
        [e1(:,:,loops),e1_pat(:,:,loops),e1_ns(:,:,loops),e1_ds(:,:,loops),e1_os(:,:,loops),Data1(loops,:),f1(loops,:)]=DW_Test(phin_path1,phin_sweep1);
        fprintf("*****************Iteration %d Sweep2*****\n",loops)
        [e2(:,:,loops),e2_pat(:,:,loops),e2_ns(:,:,loops),e2_ds(:,:,loops),e2_os(:,:,loops),Data2(loops,:),f2(loops,:)]=DW_Test(phin_path2,phin_sweep2);
    end

    %turn off Phi+, Phi-, and FFL
    if Phin, c_ma(icb_chan,0); else, c_ma(icb_chan,[0,0]); end
    if FF, c_ma(ffl_chan,0); end   
    
    disp("***********************************************************************")
    %temp grab
    if contains(fridge, "he", "IgnoreCase", true)
        temp_string="4.2K He ";
    elseif strcmpi(fridge,"adr1") || strcmpi(fridge,"adr")
        temp_string=string(round(adr_temp*1000))+"mK ";
    elseif strcmpi(fridge,"adr2")
        temp_string=string(round(adr2_temp*1000))+"mK ";
    else
        temp_string=" ";
    end

    %Add tag for DUT based in icb_chan name
    icb_char=char(string(icb_chan(1)));
    und_idx=find(icb_char == '_', 1);
    if ~isempty(und_idx), dev_string=string(icb_char(1:und_idx-1)); else, dev_string=""; end

    %Add tag to device string labeling test order
    switch test_flag
        case 1
            dev_string=string(dev_string)+ " ER ";
        case 2
            dev_string=string(dev_string)+ " VNA ";
        case 3
            dev_string=string(dev_string)+ " VNA,ER ";
        case 4
            dev_string=string(dev_string)+ " ER,VNA ";
    end

    %Plot DW of comp ER
    if test_flag~=2  
        cutoff=0.5;
        e1_cutoff=e1>cutoff;
        e2_cutoff=e2>cutoff;
        e1_cutoff_mean=mean(mean(e1_cutoff,1),3);
        e2_cutoff_mean=mean(mean(e2_cutoff,1),3);
        e1_mean=mean(mean(e1,1),3);
        e2_mean=mean(mean(e2,1),3);
        e1_pat_mean=mean(mean(e1_pat,1),3);
        e2_pat_mean=mean(mean(e2_pat,1),3);
        e1_ds_mean=mean(mean(e1_ds,1),3);
        e2_ds_mean=mean(mean(e2_ds,1),3);
        e1_os_mean=mean(mean(e1_os,1),3);
        e2_os_mean=mean(mean(e2_os,1),3);

        if (ischar(bitwise) || isstring(bitwise)) && contains(lower(string(bitwise)),"all")
            bit_visibility="on";
            pattern_visibility="on";
            cutoff_visibility="on";
            ds_os_visibility="on";
            fig_name="DW All ERs v "+Phin_string+" "+string(dev_string)+", at "+string(FFL_or_Phip)+"="+PI_string+"mA "+ temp_string+string(ntest)+"xall loopedx"+string(loop_times);
        elseif (ischar(bitwise) || isstring(bitwise)) && contains(lower(string(bitwise)),"cutoff")
            bit_visibility="off";
            pattern_visibility="off";
            cutoff_visibility="on";
            ds_os_visibility="off";
            fig_name="DW Pat Cutoff ER v "+Phin_string+" "+string(dev_string)+", at "+string(FFL_or_Phip)+"="+PI_string+"mA "+ temp_string+string(ntest)+"x"+string(bits)+"b"+" loopedx"+string(loop_times);
        elseif ((isnumeric(bitwise) && bitwise<=0) || (islogical(bitwise) && ~bitwise)) || ...
                ((ischar(bitwise) || isstring(bitwise)) && contains(lower(string(bitwise)),"pat"))
            bit_visibility="off";
            pattern_visibility="on";
            cutoff_visibility="off";
            ds_os_visibility="off";
            fig_name="DW Pat ER v "+Phin_string+" "+string(dev_string)+", at "+string(FFL_or_Phip)+"="+PI_string+"mA "+ temp_string+string(ntest)+"x"+string(bits)+"b"+" loopedx"+string(loop_times);
        else
            bit_visibility="on";
            pattern_visibility="off";
            cutoff_visibility="off";
            ds_os_visibility="off";
            fig_name="DW Bit ER v "+Phin_string+" "+string(dev_string)+", at "+string(FFL_or_Phip)+"="+PI_string+"mA "+ temp_string+string(ntest)+"x"+string(bits)+"b"+" loopedx"+string(loop_times);       
        end
        figure
        plot(phin_sweep1,e1_mean,"DisplayName", "Sweep1 Bit ER","Visible",bit_visibility,"HandleVisibility",bit_visibility)
        hold on
        plot(phin_sweep2,e2_mean,"DisplayName", "Sweep2 Bit ER","Visible",bit_visibility,"HandleVisibility",bit_visibility)
        plot(phin_sweep1,e1_pat_mean,"DisplayName", "Sweep1 Full Pat ER","Visible",pattern_visibility,"HandleVisibility",pattern_visibility)
        plot(phin_sweep2,e2_pat_mean,"DisplayName","Sweep2 Full Pat ER","Visible",pattern_visibility,"HandleVisibility",pattern_visibility)
        plot(phin_sweep1,e1_cutoff_mean,"DisplayName", "Sweep1 Cutoff Pat ER", "Visible", cutoff_visibility,"HandleVisibility",cutoff_visibility)
        plot(phin_sweep2,e2_cutoff_mean,"DisplayName", "Sweep2 Cutoff Pat ER", "Visible", cutoff_visibility,"HandleVisibility",cutoff_visibility)
        plot(phin_sweep1,e1_ds_mean+e1_os_mean,"DisplayName", "Sweep1 DS+OS BER", "Visible", ds_os_visibility,"HandleVisibility",ds_os_visibility)
        plot(phin_sweep2,e2_ds_mean+e2_os_mean,"DisplayName", "Sweep2 DS+OS BER", "Visible", ds_os_visibility,"HandleVisibility",ds_os_visibility)
        xlabel("Phi- (mA)")
        ylabel("Error Rate")
        legend
        title({string(fig_name),string(tag)})
        saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
        hold off
        close(gcf)
    end

    %Plot DW of VNA S21
    if test_flag~=1
        %Average VNA data across loop iterations while retaining the individual results.
        Data1_mean = cell(1, size(Data1, 2));
        f1_mean = cell(1, size(f1, 2));
        Data2_mean = cell(1, size(Data2, 2));
        f2_mean = cell(1, size(f2, 2));
        for col = 1:size(Data1, 2)
            Data1_mean{col} = mean(cat(3, Data1{:, col}), 3, 'omitnan');
            f1_mean{col} = mean(cat(3, f1{:, col}), 3, 'omitnan');
            Data2_mean{col} = mean(cat(3, Data2{:, col}), 3, 'omitnan');
            f2_mean{col} = mean(cat(3, f2{:, col}), 3, 'omitnan');
        end
        
        %Plot 1D DW of VNA S21 at specified frequency(ies)
        %rang=linspace(fstart,fstop,npoints);
        figure
        hold on
        for k=1:size(sample_range,1)            
            mini=sample_range(k,1);
            maxi=sample_range(k,2);
            idx1=find(rang>=mini,1,"first");
            idx2=find(rang<=maxi,1,"last");
            if isempty(idx1), idx1=1; end
            if isempty(idx2), idx2=length(rang); end
            AvgPwr1={};
            AvgPwr2={};
            for j=1:length(phin_sweep1) 
                AvgPwr1{j,k}=[mean(Data1_mean{j}(idx1:idx2,2)),   mean(atan(Data1_mean{j}(idx1:idx2,4)./Data1_mean{j}(idx1:idx2,3))),    mean(Data1_mean{j}(idx1:idx2,3)),    mean(Data1_mean{j}(idx1:idx2,4))];
                AvgPwr2{j,k}=[mean(Data2_mean{j}(idx1:idx2,2)),   mean(atan(Data2_mean{j}(idx1:idx2,4)./Data2_mean{j}(idx1:idx2,3))),    mean(Data2_mean{j}(idx1:idx2,3)),    mean(Data2_mean{j}(idx1:idx2,4))];
            end            
            p1=cell2mat(AvgPwr1(:,k));
            p2=cell2mat(AvgPwr2(:,k));
            plot(phin_sweep1, p1(:,1:4:end),"DisplayName", "dB Swp1 "+string(round(mean(sample_range(k,:))./1e9,2)))
            plot(phin_sweep2, p2(:,1:4:end),"DisplayName", "dB Swp2 "+string(round(mean(sample_range(k,:))./1e9,2)))
            plot(phin_sweep1, p1(:,2:4:end),"DisplayName", "Phase Swp1 "+string(round(mean(sample_range(k,:))./1e9,2)), "Visible", "off", "YLimInclude","off","HandleVisibility","off")       
            plot(phin_sweep2, p2(:,2:4:end),"DisplayName", "Phase Swp2 "+string(round(mean(sample_range(k,:))./1e9,2)), "Visible", "off", "YLimInclude","off","HandleVisibility","off") 
            plot(phin_sweep1, p1(:,3:4:end),"DisplayName", "Re(S21) Swp1 "+string(round(mean(sample_range(k,:))./1e9,2)), "Visible", "off", "YLimInclude","off","HandleVisibility","off") 
            plot(phin_sweep2, p2(:,3:4:end),"DisplayName", "Re(S21) Swp2 "+string(round(mean(sample_range(k,:))./1e9,2)), "Visible", "off", "YLimInclude","off","HandleVisibility","off") 
            plot(phin_sweep1, p1(:,4:4:end),"DisplayName", "Im(S21) Swp1 "+string(round(mean(sample_range(k,:))./1e9,2)), "Visible", "off", "YLimInclude","off","HandleVisibility","off") 
            plot(phin_sweep2, p2(:,4:4:end),"DisplayName", "Im(S21) Swp2 "+string(round(mean(sample_range(k,:))./1e9,2)), "Visible", "off", "YLimInclude","off","HandleVisibility","off")                 
        end
        axis tight
        hold off
        xlabel("Phi- (mA)")
        ylabel("S21")
        legend
        fig_name="DW S21 dB v "+Phin_string+" " +string(dev_string)+", at "+string(FFL_or_Phip)+"="+PI_string+"mA, "+temp_string+string(round(mean(sample_range(1,:))./1e9,2)) + "GHz"+" loopedx"+string(loop_times); 
        title({string(fig_name),string(tag)})
        saveas(gcf,string(fig_name) + " " +string(tag)+".fig")
        close(gcf)

        %Plot 2D Heatmaps of DW S21 vs Phi- & Freq if frequency span more than 1MHz, 
        if abs(fstop-fstart) >.001e9      
            div=ceil((npoints-1)/100)*10+1;
            ylabels=string(round(f1_mean{1}'./1e6));
            ylabels(  mod(1:length(f1_mean{1}),div)~=1   ) ="";
            xlabels1=string(round(phin_sweep1(:),3));
            xlabels1(  mod(1:length(phin_sweep1),5)~=1   ) ="";
            c1=cell2mat(Data1_mean(:)');
            xlabels2=string(round(phin_sweep2(:),3));
            xlabels2(  mod(1:length(phin_sweep2),5)~=1   ) ="";
            c2=cell2mat(Data2_mean(:)');
        
            figure
            h=heatmap(flip(phin_sweep1), f1_mean{1}'./1e6, flip(c1(:,2:4:end),2), "GridVisible", "off", "XDisplayLabels", flip(xlabels1),"YDisplayLabels", ylabels, "ColorLimits", [min(c1(:,2:4:end),[],"all"), max(c1(:,2:4:end),[],"all")], "Colormap",jet);
            hs = struct(h);
            xlabel("Phi- (mA)");
            ylabel("Freq (MHz)");
            ylabel(hs.Colorbar,"dB");
            fig_name="DW Freq Htmp S21 dB v "+Phin_string+" "+string(dev_string)+", at "+string(FFL_or_Phip)+"="+PI_string+"mA Swp1, "+temp_string+string(dev_string)+"loopedx"+string(loop_times); 
            title({string(fig_name),string(tag)})
            saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
            close(gcf)
            
            figure
            h=heatmap(phin_sweep2, f2_mean{1}'./1e6, c2(:,2:4:end), "GridVisible", "off", "XDisplayLabels", xlabels2,"YDisplayLabels", ylabels, "ColorLimits", [min(c2(:,2:4:end),[],"all"), max(c2(:,2:4:end),[],"all")], "Colormap",jet);
            hs = struct(h);
            xlabel("Phi- (mA)");
            ylabel("Freq (MHz)");
            ylabel(hs.Colorbar,"dB");
            fig_name="DW Freq Htmp S21 dB v "+Phin_string+" "+string(dev_string)+", at "+string(FFL_or_Phip)+"="+PI_string+"mA Swp2, "+temp_string+string(dev_string)+"loopedx"+string(loop_times);   
            title({string(fig_name),string(tag)})
            saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
            close(gcf)
         %Code for heatmaps of S21 Phase vs Phi- & Freq, dont plot
            %{
            figure 
            h=heatmap(flip(phin_sweep1), f1_mean{1}'./1e6, flip(atan(c1(:,4:4:end)./c1(:,3:4:end)),2), "GridVisible", "off", "XDisplayLabels", flip(xlabels1),"YDisplayLabels", ylabels, "ColorLimits", [-1.25,1.25], "Colormap",turbo);
            hs = struct(h);
            xlabel("Phi- (mA)");
            ylabel("Freq (MHz)");
            ylabel(hs.Colorbar,"Phase");
            fig_name="DW Freq Htmp S21 Phase v "+Phin_string+" at "+string(FFL_or_Phip)+"="+PI_string+"mA Swp1, "+temp_string+string(dev_string)+"loopedx"+string(loop_times); 
            title({string(fig_name),string(tag)})
            saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
            close(gcf)
        
            figure 
            h=heatmap(phin_sweep2, f2_mean{1}'./1e6, atan(c2(:,4:4:end)./c2(:,3:4:end)), "GridVisible", "off", "XDisplayLabels", xlabels2,"YDisplayLabels", ylabels, "ColorLimits", [-1.25,1.25], "Colormap",turbo);
            hs = struct(h);
            xlabel("Phi- (mA)");
            ylabel("Freq (MHz)");
            ylabel(hs.Colorbar,"Phase");
            fig_name="DW Freq Htmp S21 Phase v "+Phin_string+" at "+string(FFL_or_Phip)+"="+PI_string+"mA Swp2, "+temp_string+string(dev_string)+"loopedx"+string(loop_times);   
            title({string(fig_name),string(tag)})
            saveas(gcf,string(fig_name) + " " +string(tag)+".fig");
            close(gcf)
            %}
        end
    end  
   %sweep test for both phip and ffl, use for sweep up and sweep down
   function  [e_test, e_pat_test,  e_ns_test, e_ds_test, e_os_test, Data_test, f_test]= DW_Test(phin_path_arr, phin_sweep_arr)
        e_test=[];      %Temporary function variable initialization 
        e_pat_test=[];
        e_ns_test=[];
        e_ds_test=[];
        e_os_test=[];
        f_test={};
        Data_test={};
        
        if Phin 
            start= phin0;                                                           %Start location for Phi- with Phi+/FFL in zero state with 1-line Phi-    
            test_sweep=[phin_sweep_arr]';                                           %define Sweep array for actual DW test points along 1-line Phi-
        else
            start= [(phip0 + phin0)./2 , (phin0 - phip0)./2];                       %Start location for Phi- with Phi+/FFL in zero state with 2-line Phi-    
            test_sweep=[(phip_PI+phin_sweep_arr)./2; (phin_sweep_arr-phip_PI)./2]'; %define Sweep array for actual DW test points along 2-line Phi-
        end
               
        for j=1:length(test_sweep)                                              %For each test point, conduct full initialization sweeps
            c_ma(icb_chan,start);                                               %Set Phip and Phi- to start state
            if FF, c_ma(ffl_chan,ffl0); end                                     %If using FFL, set in 0 state to start

            if Phin 
                sweep=[phin_path_arr]';                                         %define sweep array for initial 1-line Phi- sweep with Phi+ (or FFL) in zero state
            else
                sweep=[(phip0+phin_path_arr)./2; (phin_path_arr-phip0)./2]';    %define sweep array for initial 2-line Phi- sweep with Phi+ (or FFL) in zero state
            end

            for i=1:length(sweep)                                               %Sweep Phi- with Phi+ or FFL in zero state to Phi- to pivot point
                c_ma(icb_chan,sweep(i,:));     
                %fprintf("Init path Meas. Phi+ = %.4f, Phi- = %.4f\n", c_ma(icb_chan(1))-c_ma(icb_chan(2)), c_ma(icb_chan(1))+c_ma(icb_chan(2))) %DO NOT RUN!!! Measuring current kicks JDPD out of state
            end
            
            if Phin 
                c_ma(icb_chan,phin_path_arr(end));                                                  %1-line Phi- at pivot point
            else
                c_ma(icb_chan,[(phip_PI+phin_path_arr(end))./2; (phin_path_arr(end)-phip_PI)./2]'); %Flip phi+ to PI state with 2-line Phi- at pivot point
            end
            if FF, c_ma(ffl_chan,ffl_PI); end                                                       %If using FFL, flip to PI state
            
            for i=1:j                                                                               %Sweep Phi- along Phi+ or FFL = PI up to test point
                c_ma(icb_chan, test_sweep(i,:));
                %fprintf("Test path Meas. Phi+ = %.4f, Phi- = %.4f\n", c_ma(icb_chan(1))-c_ma(icb_chan(2)), c_ma(icb_chan(1))+c_ma(icb_chan(2))) %DO NOT RUN!!! Measuring current kicks JDPD out of state
            end


            if Phin
                fprintf("FFL = %.4f, Phi- = %.4f\n", ffl_PI , phin_sweep_arr(j))
                %fprintf("Test Point Meas. FFL = %.4f, Phi- = %.5f\n", c_ma(ffl_chan), c_ma(icb_chan)) %DO NOT RUN!!! Measuring current kicks JDPD out of state
            else
                if FF 
                    fprintf("FFL = %.4f, Phi- = %.4f\n", ffl_PI , phin_sweep_arr(j))
                    %fprintf("Test Point Meas. FFL = %.4f, Phi- = %.4f\n", c_ma(ffl_chan), c_ma(icb_chan(1))+c_ma(icb_chan(2)))  %DO NOT RUN!!! Measuring current kicks JDPD out of state 
                elseif ~FF 
                    fprintf("Phi+ = %.4f, Phi- = %.4f\n", phip_PI , phin_sweep_arr(j))
                    %fprintf("Test Point Meas. Phi+ = %.4f, Phi- = %.4f\n", c_ma(icb_chan(1))-c_ma(icb_chan(2)), c_ma(icb_chan(1))+c_ma(icb_chan(2))) %DO NOT RUN!!! Measuring current kicks JDPD out of state
                end
            end

            switch test_flag                                                                 %Run tests with Phi- at test point with Phi+ or FFL = PI in specified order
                case 1
                    [e_test(:,j),e_pat_test(:,j),e_ns_test(:,j),e_ds_test(:,j),e_os_test(:,j)]=test_ER();
                case 2
                    [Data_test{j}, f_test{j}]=test_VNA();
                case 3
                    [Data_test{j}, f_test{j}]=test_VNA();
                    [e_test(:,j),e_pat_test(:,j),e_ns_test(:,j),e_ds_test(:,j),e_os_test(:,j)]=test_ER();
                case 4
                    [e_test(:,j),e_pat_test(:,j),e_ns_test(:,j),e_ds_test(:,j),e_os_test(:,j)]=test_ER();
                    [Data_test{j}, f_test{j}]=test_VNA();
           end
        end
        if test_flag==1
            Data_test= {nan};
            f_test= {nan};
        elseif test_flag==2
            e_test=nan;
            e_pat_test=nan;   
            e_ns_test=nan;  
            e_ds_test=nan;  
            e_os_test=nan;  
        end
    end
    function [e_bit_arg, e_pat_arg, e_ns_arg, e_ds_arg, e_os_arg]=test_ER()     
		e_pat_arg= [];
        e_bit_arg= [];
		e_ns_arg= [];
		e_ds_arg= [];
		e_os_arg= [];
        for n=1:ntest  %test pattern ntimes at each bias level
            [p, b_er, ns_er, ds_er, os_er] = pattern_test_spliterrors(pattest,edges,[],0);  %run test pattern, no patclears, assume trigger on variable "edges" (defined as 1 or rising edge), split and categorize errors
			e_pat_arg(n,1)= (1-p);
			e_bit_arg(n,1)= b_er/bits;
			e_ns_arg(n,1)=  ns_er/bits;
			e_ds_arg(n,1)=  ds_er/bits;
			e_os_arg(n,1)=  os_er/bits;
        end
		fprintf("Pat ER is %.2f and Bit ER rate is %.4f [%.2f NS %.2f DS %.2f OS]  \n", mean(e_pat_arg), mean(e_bit_arg), mean(e_ns_arg), mean(e_ds_arg), mean(e_os_arg));
    end  
    function [Data_arg, f_arg]=test_VNA()
        [f_arg, dataF_arg, dataS_arg]= v.getData;
        Data_arg=[f_arg' dataF_arg' dataS_arg(1:2:end)'  dataS_arg(2:2:end)']; 
        avgp=mean(Data_arg(idx1:idx2,2));
        fprintf("@%.2f GHz, Avg S21 = %.2f dB\n", round(mean(sample_range(1,:))./1e9,2), avgp)       
    end
end