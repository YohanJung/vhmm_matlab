function [se sp TP TN FP FN] = arcroc(True_arcs,VB_arcs,CSignBool)
% Receiver Operator Characteristics of arcs
% Compare the infered network using VBSSM with the real genetic network

% Input:
% True_arcs -- inter matrix of true network
% VB_arcs   -- inter matrix of infered network by VBSSM
% CSignBool -- 1: care the sign of weights (default)
%              0: don't care the sign of weights
% Output:
% se: sensitivity, that is, the proportion of recovered true arcs.
% sp: specificity, that is, the proportion of erroneously recovered spurious arcs.
% TP: number of true  positive
% TN: number of true  negative
% FP: number of false positive
% FN: number of false negative
%
% Rule to determine Positive/Negative and True/False

% Positive -- Presence of arc between 2 nodes in infered network
% Negative -- Abscence of arc between 2 nodes in infered network
% True     -- Both the existence and direction of arc(i-->j) in infered network
%             agree with arc(i-->j) in true network
% False    -- Existence or direction of arc(i-->j) in infered network
%             doesn't agree with arc(i-->j) in true network

%-------------------------------------------------------------------------
%         Infered   |  Real   |                Classification
%--------------------------------------------------------------------------
%   Weight [Exist?] |  Weight |[Sign Agree?]|[care sign]|| [not care sign]
%--------------------------------------------------------------------------
%        -1    [P]  |   -1    |     [T]     |   TP*     ||   TP*
%        -1    [P]  |    0    |     [F]     |   FP      ||   FP
%        -1    [P]  |    1    |     [F]     |   FP      ||   TP*
%--------------------------------------------------------------------------
%         0    [N]  |   -1    |     [F]     |   FN      ||   FN
%         0    [N]  |    0    |     [T]     |   TN*     ||   TN
%         0    [N]  |    1    |     [F]     |   FN      ||   FN
%--------------------------------------------------------------------------
%         1    [P]  |   -1    |     [F]     |   FP      ||   TP*
%         1    [P]  |    0    |     [F]     |   FP      ||   FP
%         1    [P]  |    1    |     [T]     |   TP*     ||   TP*
%--------------------------------------------------------------------------

% Juan Li, SUNY-Buffalo CSE, 10/22/2005
% juanli@cse.buffalo.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%By default, we care for the direction of arcs.
if nargin<3, CSignBool = 1;end

p = size(VB_arcs,1);

if size(True_arcs)~=size(VB_arcs)
    error('matrix sizes mismatch!');
end

TP = 0; TN = 0; FN = 0; FP = 0;
if CSignBool
    for i = 1:p
        for j = 1:p
            % Presence in both networks, AND same direction
            if True_arcs(i,j)*VB_arcs(i,j)>0
                TP = TP + 1;
            end
            %Absence in both infered and true network
            if ~VB_arcs(i,j) & ~True_arcs(i,j)
                TN = TN + 1;
            end
            % Absence in infered network, BUT presence in true network
            if ~VB_arcs(i,j) & True_arcs(i,j)
                FN = FN + 1;
            end
            %Presence in both networks, BUT different direction
            if VB_arcs(i,j) & True_arcs(i,j)*VB_arcs(i,j)<=0
                FP = FP + 1;
            end
        end
    end
else
    for i = 1:p
        for j = 1:p
            % Presence in both networks, AND same direction
            if True_arcs(i,j)& VB_arcs(i,j)
                TP = TP + 1;
            end
            %Absence in both infered and true network
            if ~VB_arcs(i,j) & ~True_arcs(i,j)
                TN = TN + 1;
            end
            % Absence in infered network, BUT presence in true network
            if ~VB_arcs(i,j) & True_arcs(i,j)
                FN = FN + 1;
            end
            % Presence in infered network, BUT absenbce in true network
            if VB_arcs(i,j) & ~True_arcs(i,j)
                FP = FP + 1;
            end
        end
    end
end

% se: proportion of recovered true arcs.
se = TP/(TP+FN);
% 1-sp:  proportion of erroneously recovered spurious arcs.
sp = TN/(TN+FP);
% complementary sp
% cmpl_sp = 1 - sp;
%
% figure,
% plot(cmpl_sp,se,'r*');
% axis([-0.1 1.1 -0.1 1.1]);
% xlabel('complementary specificity');
% ylabel('sentitivity');
%

