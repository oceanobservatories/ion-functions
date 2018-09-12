% This m-file reads raw data from Sunburst SAMI pCO2 data file and
% calculates pCO2.
% Code supplied by Reggie Spaulding at Sunburst Sensors in Feb 2018.
%
% Test data file was created from data collected during the deployment of
% sensor SN C0123 in July 2017.
%
% *************************************************************************
clear all

% 
% Sensor Specific Calibration coefficients for C0123, 2017
%
CalT=4.6539;
CalA=0.0422;
CalB=0.6761;
CalC=-1.5798;

% e's for calculating pCO2, constants provided by Sunburst based on lab
% determinations with BTB
e1=0.0043;
e2=2.136;
e3=0.2105;

File = 'testFile.txt';
fid = fopen(File);

% Read in SAMI hex data
i=1; j=1; k=1;
while 1
    s=fgetl(fid);
    if s == -1, break,end % indicates EOF
    if ~isempty(s)
        if (strcmpi('*',s(1:1))==1 && strcmpi('04',s(6:7))| strcmpi('05',s(6:7)))
            if (length(s) == 81)
                s = s(2:length(s));
                while length(s) < 80
                    s = [s,fgetl(fid)];
                end
                AA(i,:)= s;
                i=i+1;
            end
        else
        end
    end
end
fclose(fid);
[d1,d2]=size(AA);
% *************************************************************************
% Extract data from hex string
for i=1:d1
    Type(i)=hex2dec(AA(i,5:6));  % Type 4 - Measurement, Type 5 - Blank
    Time(i)=hex2dec(AA(i,7:14));  % Time
    t_SAMI(i) = ((Time(i)/86400) +(365.25*4)+(1)); 
    DarkRef(i)=hex2dec(AA(i,15:18));  % Dark Reference LED
    DarkSig(i)=hex2dec(AA(i,19:22)); % Dark Signal LED
    Ref434(i)=hex2dec(AA(i,23:26));  % 434nm Reference LED intensity
    Sig434(i)=hex2dec(AA(i,27:30));  % 434nm Signal Signal LED intensity
    Ref620(i)=hex2dec(AA(i,31:34));  % 620nm Reference LED intensity
    Sig620(i)=hex2dec(AA(i,35:38));  % 434nm Signal Signal LED intensity
    Ratio434(i)=hex2dec(AA(i,39:42)); % 434nm Ratio
    Ratio620(i)=hex2dec(AA(i,43:46)); % 620nm Ratio
    Battery(i)=hex2dec(AA(i,71:74)); % Raw Battery Value
    BattV(i) = Battery(i)*15/4096; % Battery voltage
    Therm(i)=hex2dec(AA(i,75:78)); % Thermistor raw value
end

k=1;
for i=1:d1
    BFlag(i)=nan;
    if Type(i)==5 % Blank measurement
        A434Bka(k)=Ratio434(i)/16384;
        A620Bka(k)=Ratio620(i)/16384;
        blanktime(k)=((Time(i))/(60*60*24)); % Time of blank
        datetime1(i)=(Time(i)/(60*60*24));   % Time total 
        Tc_SAMI(i)=(1./((0.0010183)+(0.000241.*log((Therm(i)./...
            (4096-Therm(i))).*17400))+((0.00000015.*log((Therm(i)./...
            (4096-Therm(i))).*17400).^3))))-273.15; % Temp
        k=k+1;
    end
end

% PCO2 Calculations
i2=1;
for i=1:d1
    if Type(i)==5
        A434Bk=A434Bka(i2);
        A620Bk=A620Bka(i2);
        calcCO2_1(i)=0;
        i2=i2+1;
    else if Type(i)==4
            k434(i)=A434Bk;
            k620(i)=A620Bk;
            A434(:,i)=-log10((Ratio434(i)./16384)/k434(i));
            A620(:,i)=-log10((Ratio620(i)./16384)/k620(i));
            Ratio(:,i)=(A620(:,i))/(A434(:,i));
            datetime1(i)=(Time(i)./(60*60*24));
            Tc_SAMI(i)=(1./((0.0010183)+(0.000241.*log((Therm(i)./...
                (4096-Therm(i))).*17400))+((0.00000015.*log((Therm(i)./...
                (4096-Therm(i))).*17400).^3))))-273.15; % Temp
            RCO2(:,i) = -log10((Ratio(:,i)-e1)./(e2-e3.*Ratio(:,i)));
            Tcor_RCO2i(:,i) = RCO2(:,i)+(0.008.*(Tc_SAMI(i)-CalT));
            Tcoeff(:,i) = (0.0075778)-(0.0012389.*Tcor_RCO2i(:,i))-...
                (0.00048757.*Tcor_RCO2i(:,i).^2);
            Tcor_RCO2(:,i) = RCO2(:,i)+(Tcoeff(:,i)).*(Tc_SAMI(i)-CalT);
            calcCO2_1(:,i) = (10.^(((-1*CalB)+((CalB^2)-(4*CalA.*...
                (CalC-Tcor_RCO2(:,i)))).^0.5)./(2*CalA)));
        end
    end
end
            
% *************************** Write output file ***************************
outfile_CO2=strcat(File(1:end-4),'_out.xls'); 
col_header={'Time','pCO2','Temp','BattV','DarkRef','DarkSig',...
    'Ratio434','Ref434','Sig434','Ratio620','Ref620','Sig620'};
datafinal_CO2=[t_SAMI;calcCO2_1;Tc_SAMI;BattV;DarkRef;DarkSig...
    ;Ratio434;Ref434;Sig434;Ratio620;Ref620;Sig620];
fid = fopen(outfile_CO2,'w');
fprintf(fid,['Time\t pCO2\t Temp\t BattV\t DarkRef\t DarkSig\t'...
    'Ratio434\t Ref434\t Sig434\t Ratio620\t Ref620\t Sig620\r\n']);
fmt = ['%8.4f \t %8.4f \t %8.4f \t %8.2f \t %8f \t %8f \t %8f \t'...
    '%8f \t %8f \t %8f \t %8f \t %8f \r\n'];
fprintf(fid,fmt,datafinal_CO2);
fclose(fid);
