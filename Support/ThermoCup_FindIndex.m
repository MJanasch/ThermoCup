%% FindIndex Function
% Modified from Avlant Nilsson
% Markus Janasch, Ph.D. Student
% Microbial Metabolic Engineering Group
% KTH School of Biotechnology, Science For Life Laboratory, Stockholm
% Created: 2016/06/21, Last modified: 2016/06-21

function [n] = MJanasch_FindIndex(haystack, needle)
    n=find(ismember(haystack,needle));
end
