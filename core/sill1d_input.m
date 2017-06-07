function [ROCK, SILL, FLUID, WellData] = sill1d_input(Filename)
% 1D Thermal cooling of intrusive
% COPYRIGHT(C) 2015, GeoModelling Solutions GmbH
%
% READ INPUT FILE
%
%
% Original author:    Karthik Iyer
% Last committed:     $Revision: 0 $
% Last changed by:    $Author: karthik $
% Last changed date:  $Date: 2012-03-07 16:05:55 +0100 (Wed, 07 Mar 2012) $
%--------------------------------------------------------------------------

%% Read Data
if ~exist(Filename, 'file')
   error('Input excel file does not exist in directory');
end

[FL_NUM,~,~]            = xlsread(Filename, 'Fluid');
[WellData.Toc,~,~]      = xlsread(Filename, 'TOC Data');
[LITH_NUM,LITH_NAME,~]  = xlsread(Filename, 'Lithology');
[SILL_NUM,~,~]          = xlsread(Filename, 'Sills');
[WellData.VR,~,~]       = xlsread(Filename, 'Vitrinite Data');
[WellData.T,~,~]        = xlsread(Filename, 'Temperature Data');
[ERO_NUM,ERO_NAME,~]    = xlsread(Filename, 'Erosion');

FLUID.rho       = FL_NUM(1,1);
FLUID.cp        = FL_NUM(1,2);

if ~isempty(ERO_NUM) 
    SORT_NUM        = [LITH_NUM;ERO_NUM(:,1:10)];
    SORT_NAME       = [LITH_NAME(2:end,1);ERO_NAME(2:end,1)];
else
    SORT_NUM        = LITH_NUM;
    SORT_NAME       = LITH_NAME(2:end,1);
end

% Sort Lithology according to Depth
[~,b]           = sort(SORT_NUM(:,2));
SORT_NUM        = SORT_NUM(b,:);

ROCK.Name       = SORT_NAME(b,:);
ROCK.num        = size(SORT_NUM,1);
ROCK.top        = SORT_NUM(1,1);
ROCK.bot        = SORT_NUM(end,1);
ROCK.Tops       = SORT_NUM(:,1);
ROCK.Ages       = SORT_NUM(:,2);
ROCK.Rho        = SORT_NUM(:,3);
ROCK.Cp         = SORT_NUM(:,4);
ROCK.Phi        = SORT_NUM(:,5);
ROCK.K          = SORT_NUM(:,6);
ROCK.Toc        = SORT_NUM(:,7)./100;
ROCK.Lm         = SORT_NUM(:,8).*1e3;
ROCK.Ld         = SORT_NUM(:,9).*1e3;
ROCK.Carb       = SORT_NUM(:,10);

ROCK.Ero_t      = NaN*ones(size(ROCK.Ld));
ROCK.Ero_thick  = NaN*ones(size(ROCK.Ld));

[~,b]           = sort(SILL_NUM(:,1));
SILL_NUM        = SILL_NUM(b,:);

SILL.num        = size(SILL_NUM,1);
SILL.Tops       = SILL_NUM(:,1);
SILL.Thick      = SILL_NUM(:,2);
SILL.E_time     = SILL_NUM(:,3);
SILL.E_temp     = SILL_NUM(:,4);
SILL.Rhom       = SILL_NUM(:,5);
SILL.Cpm        = SILL_NUM(:,6);
SILL.Rhos       = SILL_NUM(:,7);
SILL.Cps        = SILL_NUM(:,8);
SILL.K          = SILL_NUM(:,9);
SILL.Sol        = SILL_NUM(:,10);
SILL.Liq        = SILL_NUM(:,11);
SILL.Lc         = SILL_NUM(:,12).*1e3;

% Inserts Eroded Layers at correct depths in reconstructed lithological
% column and move layer below by required amounts
if ~isempty(ERO_NUM)
    Ero_age       = sort(ERO_NUM(:,2), 'descend');
    id_all        = NaN*ones(length(Ero_age),1);
    for i = 1:length(Ero_age)
        id                  = find(ROCK.Ages==Ero_age(i));
        q                   = find(ROCK.Tops==ROCK.Tops(id));
        q(q==id)            = [];
        if ROCK.Ages(q)<ROCK.Ages(id)
            error('sill1d_input:eroded_ages', 'Error.\nEroded Layer Age cannot be younger than the Top Age of the Lithology at the same interface')
        end        
        id_all(length(Ero_age)+1-i) = id;
        ib                  = find(ERO_NUM(:,2)==Ero_age(i));
        ROCK.Ero_t(id)      = ERO_NUM(ib,11);
        ROCK.Ero_thick(id)  = ERO_NUM(ib,12);
        if id<length(ROCK.Tops) && id~=1
            ROCK.Tops(id+1:end)                 = ROCK.Tops(id+1:end) + ROCK.Ero_thick(id);
            SILL.Tops(SILL.Tops>ROCK.Tops(id))  = SILL.Tops(SILL.Tops>ROCK.Tops(id)) + ROCK.Ero_thick(id);
        elseif id==length(ROCK.Tops)
            SILL.Tops(SILL.Tops>ROCK.Tops(id))  = SILL.Tops(SILL.Tops>ROCK.Tops(id)) + ROCK.Ero_thick(id);
        elseif id == 1
            ROCK.Tops(id+1:end)                 = ROCK.Tops(id+1:end) + ROCK.Ero_thick(id);
            SILL.Tops(SILL.Tops>ROCK.Tops(id))  = SILL.Tops(SILL.Tops>ROCK.Tops(id)) + ROCK.Ero_thick(id);            
        end
    end
    ROCK.Ero_tops                       = ROCK.Tops(id_all);
end


%% Error checking
% Check if Sill Tops coincide with Rock Tops
% Impossible situation
if any(ismember(SILL.Tops, ROCK.Tops))
    error('sill1d_input:sill_top_error', 'Error.\nSill tops are not allowed to coincide with rock tops.');
end

% Check if any rock tops are contained inside sills
for i = 1:length(SILL.Tops)    
    if any( ROCK.Tops>SILL.Tops(i) & (SILL.Tops(i)+SILL.Thick(i))>ROCK.Tops ) 
        error('sill1d_input:rock_top_inside_sill', 'Error.\nRock top inside sill.');
    end
end

% Checks if sills cross each other and returns error
[Stops, Sort_ind] = sort(SILL.Tops);
Sbots = Stops + SILL.Thick(Sort_ind);
if any(Sbots(1:end-1)-Stops(2:end)>0)
    error('sill1d_input:sill_crossing', 'Error.\nSills cannot cross each other.');
end

% Checks if deeper rocks have earlier ages
if any(ROCK.Ages(2:end)-ROCK.Ages(1:end-1)<0)
    error('sill1d_input:rock_ages', 'Error.\nDeeper rocks have younger ages than those above.');
end

% Checks if time of erosion is later than deposition time
if any((ROCK.Ages-ROCK.Ero_t)<=0)
    error('sill1d_input:eroded_time', 'Error.\nErosion of layer occurs before or at the same time as deposition.');
end

% Checks if older rocks are eroded after younger rocks
ind_ero = find(~isnan(ROCK.Ero_t));
for i = 1:length(ind_ero)-1
    if ROCK.Ero_t(ind_ero(i+1))-ROCK.Ero_t(ind_ero(i))>0 && ROCK.Ero_t(ind_ero(i+1))<=ROCK.Ages(ind_ero(i))
         error('sill1d_input:eroded_time', 'Error.\nErosion of layer occurs after erosion of older layer.');
    end
end

% Checks if all depositional, sill emplacement and erosion events are unique
a = length([ROCK.Ages; unique(SILL.E_time); ROCK.Ero_t]);
b = length(unique([ROCK.Ages; SILL.E_time; ROCK.Ero_t]));
if a~=b
    error('sill1d_input:unique_events', 'Error.\nTop Ages, Erosion and Sill Emplacement times are non-unique.');
end

if size(WellData.T,1)<2
    error('sill1d_input:temperature_data', 'Error.\nTemperature data must contain at least two data points with the first data point describing surface temperature.');
end