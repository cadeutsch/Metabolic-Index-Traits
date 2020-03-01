function [ Par ] = taxon_categ( Par, phynames )

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Phylum names

phyall=Par.Phylum;
phyabb=phyall;
phyunq=unique(phyall);

% All Phyla/Categories
phynumAll=zeros(1,length(phyall));
for i = 1:length(phyunq)
    I=strcmp(phyunq(i),phyall);
    phynumAll(I)=i;
end

% Abbreviated Categories
phynumAbb=zeros(1,length(phyall));
for i=1:(length(phynames)-1)
    I=strcmp(phynames{i},phyall);
    phynumAbb(I)=i;
end
I=phynumAbb==0;
phynumAbb(I)=length(phynames);
phyabb(I)=phynames(end);

% bundle into structure
Par.categ=phyabb;
Par.categ_num=phynumAbb';
Par.phyla=phyall;
Par.phyla_num=phynumAll';
Par.phyla_unq=phyunq;



end

