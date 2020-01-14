% [xcut, fs]= cutfile_generalized(name, heure, duree)
%
% Author: Lea Bouffaut, PhD
% Naval Academy Research Institute, Brest, France.
% Date: 10-30-2019
%
% ENTREES
%     name  =  nom + emplacement relatif du fichier 
%              ex. MONDOSSIER/YV.RR48.00.BDH.M.2012.318.235940.SAC.wav
%     heure =  heure de debut de la selection (heure.min)
%              ex. 00.34  21.00  10.10   08.55
%     duree =  duree souhaitee de la selection EN MINUTES
%
% SORTIES
%     xcut = vecteur du signal extrait
%     fs   = fr?quence d'echantillonnage
%
%
% EXEMPLE PROGRAMMATION
% nom   = 'MONDOSSIER/YV.RR48.00.BDH.M.2012.318.235940.SAC.wav';
% h = 11.30;
% duree = 60 ; %(min)
% acceleration = 10 ; % (*combien on souhaite accelerer)
%
% [x, fs] = cutfile_generalized(nom, h, duree) ;
%
% soundsc(x,fs*acceleration)

function [xcut, fs,N_heure_deb,N_heure_fin]= cutfile_generalized(name, heure, duree)

% vers = strcmp(version('-release'),'2016a');
% vers2 = strcmp(version('-release'),'2015b');
% if or(vers,vers2) == 1 
    [x,fs] = audioread(name);
% else
%     [x,fs] = wavread(name);
% end

x = x'; M = length(x);
% T_24h = (M-1)/fs; % Duree totale de l'enregistrement en secondes

% Coupure
T_heure = floor(heure)*3600 + (heure-floor(heure))*60*100; % D?but de l'observation en secondes
N_heure_deb = round(T_heure * fs +1 ); % Debut de l'observation en nombre d'echantillons
T_cut = duree*60; % Temps d'observation en secondes 
N_cut = T_cut *fs +1; % Temps d'observation en nombre d'echantillons
N_heure_fin = N_heure_deb + N_cut; % Fin de l'observation en nombre d'echantillons

xcut = x(N_heure_deb :N_heure_fin);
