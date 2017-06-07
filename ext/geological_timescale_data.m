function [Period_age, Period_name, Period_rgb, Stage_age, Stage_name, Stage_rgb] = geological_timescale_data()
% geological_timescale_data
% 
% Returns the data (age, name, color) of the geological periods and stages.
% Data obtained from stratigraphy.org, tscreator, and 
% http://blog.effjot.net/en/2010/11/stylesheet-stratigraphie-jetzt-komplett
%
% Note 1:
% Stages only work back to Cambrium!
%
% There is no real reason to store the color as cmyk. Should switch to rgb.
% This was the case in an original version that, however, did not produce
% the right colors due to a mistake in the corresponding pdf chart (?).
%
% Original author:    Schmid
% Last committed:     $Revision: 243 $
% Last changed by:    $Author: schmid $
% Last changed date:  $Date: 2012-03-07 16:05:55 +0100 (Wed, 07 Mar 2012) $
%--------------------------------------------------------------------------

%STRATIGRAPHY.ORG DATA
Stage_age   = [0 0.0118 0.126 0.781 1.806 2.588 3.600 5.332 7.246 11.608 13.82 15.97 20.43 23.03 28.8 33.9 37.2 40.4 48.6 55.8 58.7 61.7 65.5 70.6 83.5 85.8 89.3 93.5 99.6 112.0 125.0 130.0 136.4 140.2 145.5...
    150.8 155.7 161.2 164.7 167.7 171.6 175.6 183.0 189.6 196.5 199.6 203.6 216.5 228.0 237.0 245.0 249.7 251.0 253.8 260.4 265.8 268.0 270.6 275.6 284.4 294.6 299.0 303.9 306.5 311.7 318.1 326.4 345.3 359.2...
    374.5 385.3 391.8 397.5 407.0 411.1 416.0 418.7 421.3 422.9 426.2 428.2 436 439 443.7 445.6 455.8 460.9 468.1 471.8 478.6 488.3];

Stage_name  = {'Holocene', 'U. Pleistocene', 'M. Pleistocene', 'L. Pleistocene', 'Gelasian', 'Piacenzian', 'Zanclean', 'Messinian', 'Tortonian', 'Serravallian', 'Langhian', 'Burdigalian', 'Aquitanian', 'Chattian', 'Rupelian', 'Priabonian', 'Bartonian', 'Lutetian', 'Ypresian', 'Thanetian', 'Selandian', 'Danian', 'Maastrichtian', 'Campanian', 'Santonian', 'Coniacian', 'Turonian', 'Cenomanian', 'Albian', 'Aptian', 'Barremian', 'Hauterivian', 'Valanginian', 'Berriasian',...
    'Thitonian', 'Kimmeridgian', 'Oxfordian', 'Callovian', 'Bathonian', 'Bajocian', 'Aalenian', 'Toarcian', 'Pliensbachian', 'Sinemurgian', 'Hettangian', 'Rhaetian', 'Norian', 'Carnian', 'Ladinian', 'Anisian', 'Olenekian', 'Induan', 'Changhsingian', 'Wuchiapingian', 'Captanian', 'Wordian', 'Roadian', 'Kungurian', 'Artinskian', 'Sakmarian', 'Asselian', 'Gzhelian', 'Kasimovian', 'Moscovian', 'Bashkirian', 'Serpukhovian', 'Visean', 'Tournaisian',...
    'Famennian', 'Frasnian', 'Givetian', 'Eifelian', 'Emsian', 'Pragian', 'Lochkovian', 'Pridoli', 'Ludfordian', 'Gorstian', 'Homerian', 'Sheinwoodian', 'Telychian', 'Aeronian', 'Rhuddanian', 'Hirnantian', 'Katian', 'Sandbian', 'Darriwilian', 'Stage 3', 'Floian', 'Tremadocian'};

Stage_cmyk = [0 5 5;     0 5 15;        0 5 20;         0 5 25;         0 0 20;   0 0 25;     0 0 30;    0 0 55;    0 0 60;   0 0 65;        0 0 70;   0 0 75;      0 0 80;     0 10 30;  0 15 35;  0 20 30;    0 25 35;   0 30 40;  0 35 45;  0 25 50;   0 25 55; 0 30 55;  5 0 45;        10 0 50;   15 0 55;   20 0 60;   25 0 65;  30 0 70;  20 0 40; 25 0 45; 30 0 50;    35 0 55;    40 0 60;      45 0 65;...
    15 0 0; 20 0 0; 25 0 0;25 0 5; 30 0 5; 35 0 5; 40 0 5; 40 5 0; 50 5 0; 60 5 0; 70 5 0; 10 25 0; 15 30 0; 20 35 0; 20 45 0; 25 50 0; 30 65 0; 35 70 0; 0 25 20; 0 30 25; 0 40 35; 0 45 40; 0 50 45; 10 45 40; 10 50 45; 10 55 50; 10 60 55; 20 10 15; 25 10 15; 30 10 20; 40 10 20; 25 15 55; 35 15 55; 45 15 55;...
    5 5 20; 5 5 30; 5 10 45; 5 15 50; 10 15 50; 10 20 55; 10 25 60; 10 0 10; 15 0 10; 20 0 10; 20 0 15; 25 0 20; 25 0 15; 30 0 20; 35 0 25; 35 0 30; 40 0 35; 45 0 40; 55 0 35; 60 0 40; 75 0 45; 80 0 50];

Period_name = {'Quaternary', 'Neogene', 'Paleogene', 'Cretaceous', 'Jurassic', 'Triassic', 'Permian', 'Carboniferous', 'Devonian', 'Silurian', 'Ordovician', 'Cambrian', 'Proterozoic', 'Archean'};
Period_age  = [0  1.806        23.03      65.5        145.5         199.6       251.0       299.0      359.2            416.0       443.7       488.3         542.0       2500           4600];
Period_cmyk = [  0 0 50;       0 10 100;  0 40 60;    50 0 75;      80 0 5;     50 80 0;    5 75 75;   60 15 30;        20 40 75;   30 0 25;    100 0 60;     50 20 65;   0 80 35;       0 100 0];

% CMYK -> RGB CONVERSION
Stage_rgb   = zeros(size(Stage_cmyk));
for i=1:size(Stage_cmyk, 1)
    Stage_rgb(i,:)  = ([1-Stage_cmyk(i,1:3)./100]);
end

Period_rgb  = zeros(size(Period_cmyk));
for i=1:size(Period_cmyk, 1)
    Period_rgb(i,:)  = ([1-Period_cmyk(i,1:3)./100]);
end