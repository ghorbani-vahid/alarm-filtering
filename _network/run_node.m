function est = run_node(model,meas,truth, k,tt,tt_lmb_update_pre, nn,iter_num, estt,id, save,x_off, y_off, method)
rng(42);
% This is the MATLAB code for ALARM filter based on...
% the Labeled Multi-Bernoulli filter proposed in
% S. Reuter, B.-T. Vo, B.-N. Vo, and K. Dietmayer, "The labelled multi-Bernoulli filter," IEEE Trans. Signal Processing, Vol. 62, No. 12, pp. 3246-3260, 2014
% http://ba-ngu.vo-au.com/vo/RVVD_LMB_TSP14.pdf
% which propagates an LMB approximation of the GLMB update proposed in
% B.-T. Vo, and B.-N. Vo, "Labeled Random Finite Sets and Multi-Object Conjugate Priors," IEEE Trans. Signal Processing, Vol. 61, No. 13, pp. 3460-3475, 2013.
% http://ba-ngu.vo-au.com/vo/VV_Conjugate_TSP13.pdf
% and
% B.-N. Vo, B.-T. Vo, and D. Phung, "Labeled Random Finite Sets and the Bayes Multi-Target Tracking Filter," IEEE Trans. Signal Processing, Vol. 62, No. 24, pp. 6554-6567, 2014
% http://ba-ngu.vo-au.com/vo/VVP_GLMB_TSP14.pdf
% using an efficient implementation of the GLMB filter proposed in
% B.-T. Vo, and B.-N. Vo, "An Efficient Implementation of the Generalized Labeled Multi-Bernoulli Filter," IEEE Trans. Signal Processing, Vol. 65, No. 8, pp. 1975-1987, 2017.
% http://ba-ngu.vo-au.com/vo/VVH_FastGLMB_TSP17.pdf
%
%
% Note 1: no dynamic grouping or adaptive birth is implemented in this code, only the standard filter with static birth is given
% Note 2: the simple example used here is the same as in the CB-MeMBer filter code for a quick demonstration and comparison purposes
% Note 3: more difficult scenarios require more components/hypotheses (thus exec time)
% ---BibTeX entry
% @ARTICLE{LMB,
% author={S. Reuter and B.-T. Vo and B.-N. Vo and K. Dietmayer},
% journal={IEEE Transactions on Signal Processing},
% title={The Labeled Multi-Bernoulli Filter},
% year={2014},
% month={Jun}
% volume={62},
% number={12},
% pages={3246-3260}}
%
% @ARTICLE{GLMB1,
% author={B.-T. Vo and B.-N. Vo
% journal={IEEE Transactions on Signal Processing},
% title={Labeled Random Finite Sets and Multi-Object Conjugate Priors},
% year={2013},
% month={Jul}
% volume={61},
% number={13},
% pages={3460-3475}}
%
% @ARTICLE{GLMB2,
% author={B.-T. Vo and B.-N. Vo and D. Phung},
% journal={IEEE Transactions on Signal Processing},
% title={Labeled Random Finite Sets and the Bayes Multi-Target Tracking Filter},
% year={2014},
% month={Dec}
% volume={62},
% number={24},
% pages={6554-6567}}
%
% @ARTICLE{GLMB3,
% author={B.-N. Vo and B.-T. Vo and H. Hung},
% journal={IEEE Transactions on Signal Processing},
% title={An Efficient Implementation of the Generalized Labeled Multi-Bernoulli Filter},
% year={2017},
% month={Apr}
% volume={65},
% number={8},
% pages={1975-1987}}
%---

%=== Setup

%output variables
est.X= cell(meas.K,1);
est.N= zeros(meas.K,1);
est.L= cell(meas.K,1);



%filter parameters
filter.T_max= 100;                  %maximum number of tracks
filter.track_threshold= 1e-3;       %threshold to prune tracks
filter.H_upd= 1000;                  %requested number of updated components/hypotheses (for GLMB update)

filter.npt= 1000;                   %number of particles per track
filter.nth= 100;                    %threshold on effective number of particles before resampling (not used here, resampling is forced at every step, otherwise number of particles per track grows)

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

filter.P_G= 0.9999999;                           %gate size in percentage
filter.gate_ngrid= 10;                           %grid size (on each dimension) for digital gating based on histogram of particle weights
filter.gate_flag= 1;                             %gating on or off 1/0

est.filter= filter;

%=== Filtering

%initial prior
if k==1
tt_lmb_update= cell(0,1);      %track table for LMB (cell array of structs for individual tracks)
end

if k>1
tt_lmb_update=tt_lmb_update_pre;
end


%recursive filtering
%for k=1:meas.K
  

    %joint predict and update, results in GLMB, convert to LMB
    glmb_update= jointlmbpredictupdate(tt_lmb_update,model,filter,meas,k, tt, nn,x_off, y_off);              T_predict= length(tt_lmb_update)+model.T_birth;
    tt_lmb_update= glmb2lmb(glmb_update);                                               T_posterior= length(tt_lmb_update);
       
    %pruning, truncation and track cleanup
    tt_lmb_update= clean_lmb(tt_lmb_update,filter);                                     T_clean= length(tt_lmb_update);
    
    %state estimation
    [estt.X{k, 1}, estt.N(k, 1), estt.L{k, 1}]= extract_estimates(tt_lmb_update,model);
    est=estt;
    %display diagnostics
    display_diaginfo(tt_lmb_update,k,est,filter,T_predict,T_posterior,T_clean);

    est.tt_lmb_update=tt_lmb_update; % Vahid
%end

if k==iter_num %last iteration
    if save=="true"
        attack.scenario="none";
        attack.intensity=0;
        truth=gen_truth(model, iter_num, attack);
        plot_results(model, truth, meas, estt,id, method);
    end
end    

end



function glmb_nextupdate= jointlmbpredictupdate(tt_lmb_update,model,filter,meas,k,tt,nn,x_off, y_off)
%---generate birth tracks
tt_birth= cell(length(model.r_birth),1);                                           %initialize cell array
for tabbidx=1:length(model.r_birth)
    tt_birth{tabbidx}.r= model.r_birth(tabbidx);                                   %birth prob for birth track
    tt_birth{tabbidx}.x= gen_gms(model.w_birth{tabbidx},model.m_birth{tabbidx},model.P_birth{tabbidx},filter.npt);  %samples for birth track
    tt_birth{tabbidx}.w= ones(filter.npt,1)/filter.npt;                                                             %weights of samples for birth track
    tt_birth{tabbidx}.l= [k;tabbidx];                                              %track label
end

%---generate surviving tracks
tt_survive= cell(length(tt_lmb_update),1);                                                                              %initialize cell array
for tabsidx=1:length(tt_lmb_update)
    wtemp_predict= compute_pS(model,tt_lmb_update{tabsidx}.x).*tt_lmb_update{tabsidx}.w(:); xtemp_predict= gen_newstate_fn(model,tt_lmb_update{tabsidx}.x,'noise');      %particle prediction
    tt_survive{tabsidx}.r= sum(wtemp_predict)*tt_lmb_update{tabsidx}.r;                                                          %predicted existence probability for surviving track
    tt_survive{tabsidx}.x= xtemp_predict;                                                                                                                                   %samples for surviving track
    tt_survive{tabsidx}.w= wtemp_predict/sum(wtemp_predict);  
    tt_survive{tabsidx}.l= tt_lmb_update{tabsidx}.l;                                                                    %track label
end

%create predicted tracks - concatenation of birth and survival
tt_predict = cat(1, tt_survive, tt_birth); %tt_predict= cat(1,tt_birth,tt_survive);                                                                                %copy track table back to GLMB struct
%tt_predict= cat(1,tt_birth,tt_survive); 
%gating by tracks
if filter.gate_flag
    for tabidx=1:length(tt_predict)
        tt_predict{tabidx}.gatemeas= gate_meas_smc_idx(meas.Z{k},filter.P_G,model,tt_predict{tabidx}.x,tt_predict{tabidx}.w,filter.gate_ngrid);
    end
else
    for tabidx=1:length(tt_predict)
        tt_predict{tabidx}.gatemeas= 1:size(meas.Z{k},2);
    end
end

%precalculation loop for average survival/death probabilities
avps= zeros(length(tt_predict),1);
for tabidx=1:length(tt_predict)
    avps(tabidx)= tt_predict{tabidx}.r;
end
avqs= 1-avps;

%precalculation loop for average detection/missed probabilities
avpd= zeros(length(tt_predict),1);
for tabidx=1:length(tt_predict)
    avpd(tabidx)= tt_predict{tabidx}.w(:)'*compute_pD(model,tt_predict{tabidx}.x)+eps(0);
end
avqd= 1-avpd;

%create updated tracks (single target Bayes update)
m= size(meas.Z{k},2);                                   %number of measurements
tt_update= cell((1+m)*length(tt_predict),1);            %initialize cell array
%missed detection tracks (legacy tracks)
for tabidx= 1:length(tt_predict)
    tt_update{tabidx}= tt_predict{tabidx};              %same track table
end
%measurement updated tracks (all pairs)
allcostm= zeros(length(tt_predict),m);

%Liklihood modification (averaging over recieved updated states from
%neighboring nodes to eliminate the impact of wrong measurements. should be
%fine when all birth tracks are the same which is the default.

num_tables = numel(tt);   %tt is the neighboring estimates

%%% comment for S-ALARM
for tabsidx=1:length(tt_predict)% apply chapman kologrov on neighboring nodes paricles of corrspondinf tracks (tabsidx) which are similar due to the label definition  
    for index=1:num_tables 
    if tabsidx <= numel(tt(index).tt_lmb_update)
        xtemp_predict= gen_newstate_fn(model,tt(index).tt_lmb_update{tabsidx}.x,'noise');
        tt(index).tt_lmb_update{tabsidx}.x=xtemp_predict;
    end
    end
end 

for tabidx= 1:length(tt_predict) 
    
    for emm= tt_predict{tabidx}.gatemeas
            stoidx= length(tt_predict)*emm + tabidx; %index of predicted track i updated with measurement j is  (number_predicted_tracks*j + i)
            new_likelihood=0;
            or_likelihood= compute_likelihood(model, meas.Z{k}(:, emm), tt_predict{tabidx}.x,0, 0);
            label=tt_predict{tabidx}.l;
            
            for index=1:num_tables
        
                exists = any(cellfun(@(x) isequal(label, x.l), tt(index).tt_lmb_update));


                if exists ==1
               
                   new_likelihood= new_likelihood+ or_likelihood;
                end

            end

                
                
            
                           
                avr_likelihood = or_likelihood; 
                
                if new_likelihood > 0
                    avr_likelihood = (or_likelihood + new_likelihood) / (num_tables + 1);
                else
                    avr_likelihood = 0 * or_likelihood; 
                end

              

                
               

                if num_tables==0 %single node (ST)
                    avr_likelihood =or_likelihood;
                end

         
            w_temp= avps(tabidx)*compute_pD(model,tt_predict{tabidx}.x).*tt_predict{tabidx}.w(:).*  avr_likelihood'; x_temp= tt_predict{tabidx}.x;  %weight update for this track and this measuremnent      
            tt_update{stoidx}.x= x_temp;                                                        %particles for updated track
            tt_update{stoidx}.w= w_temp/sum(w_temp);                                            %weights of partcles for updated track
            tt_update{stoidx}.l = tt_predict{tabidx}.l;                                                                                     %track label
            allcostm(tabidx,emm)= sum(w_temp);                                                                                              %predictive likelihood
    end
end

glmb_nextupdate.tt= tt_update;                                                                                                              %copy track table back to GLMB struct
%joint cost matrix
jointcostm= [diag(avqs) ...
             diag(avps.*avqd) ...
             allcostm/(model.lambda_c*model.pdf_c)];
%gated measurement index matrix
gatemeasidxs= zeros(length(tt_predict),m);
for tabidx= 1:length(tt_predict)
    gatemeasidxs(tabidx,1:length(tt_predict{tabidx}.gatemeas))= tt_predict{tabidx}.gatemeas;
end
gatemeasindc= gatemeasidxs>0;
         

%component updates

    %calculate best updated hypotheses/components
    cpreds= length(tt_predict);
    nbirths= model.T_birth;
    nexists= length(tt_lmb_update);
    ntracks= nbirths + nexists;
    tindices= [(1:nbirths) nbirths+(1:nexists)];                                                                                          %indices of all births and existing tracks  for current component
    lselmask= false(length(tt_predict),m); lselmask(tindices,:)= gatemeasindc(tindices,:);                                              %logical selection mask to index gating matrices
    mindices= unique_faster(gatemeasidxs(lselmask));                                                                                    %union indices of gated measurements for corresponding tracks
    costm= jointcostm(tindices,[tindices cpreds+tindices 2*cpreds+mindices]);                                                           %cost matrix - [no_birth/is_death | born/survived+missed | born/survived+detected]
    neglogcostm= -log(costm);                                                                                                           %negative log cost
    [uasses,nlcost]= gibbswrap_jointpredupdt_custom(neglogcostm,round(filter.H_upd));                                                   %murty's algo/gibbs sampling to calculate m-best assignment hypotheses/components
    uasses(uasses<=ntracks)= -inf;                                                                                                      %set not born/track deaths to -inf assignment
    uasses(uasses>ntracks & uasses<= 2*ntracks)= 0;                                                                                     %set survived+missed to 0 assignment
    uasses(uasses>2*ntracks)= uasses(uasses>2*ntracks)-2*ntracks;                                                                       %set survived+detected to assignment of measurement index from 1:|Z|    
    uasses(uasses>0)= mindices(uasses(uasses>0));                                                                                       %restore original indices of gated measurements
    
    %generate corrresponding jointly predicted/updated hypotheses/components
    for hidx=1:length(nlcost)
        update_hypcmp_tmp= uasses(hidx,:)'; 
        update_hypcmp_idx= cpreds.*update_hypcmp_tmp+[(1:nbirths)'; nbirths+(1:nexists)'];
        glmb_nextupdate.w(hidx)= -model.lambda_c+m*log(model.lambda_c*model.pdf_c)-nlcost(hidx);                                             %hypothesis/component weight
        glmb_nextupdate.I{hidx}= update_hypcmp_idx(update_hypcmp_idx>0);                                                                                              %hypothesis/component tracks (via indices to track table)
        glmb_nextupdate.n(hidx)= sum(update_hypcmp_idx>0);                                                                                                            %hypothesis/component cardinality
    end

glmb_nextupdate.w= exp(glmb_nextupdate.w-logsumexp(glmb_nextupdate.w));                                                                                                                 %normalize weights

%extract cardinality distribution
for card=0:max(glmb_nextupdate.n)
    glmb_nextupdate.cdn(card+1)= sum(glmb_nextupdate.w(glmb_nextupdate.n==card));                                                                                                       %extract probability of n targets
end

%remove duplicate entries and clean track table
glmb_nextupdate= clean_update(clean_predict(glmb_nextupdate),filter);
end



function glmb_temp= clean_predict(glmb_raw)
%hash label sets, find unique ones, merge all duplicates
for hidx= 1:length(glmb_raw.w)
    glmb_raw.hash{hidx}= sprintf('%i*',sort(glmb_raw.I{hidx}(:)'));
end

[cu,~,ic]= unique(glmb_raw.hash);

glmb_temp.tt= glmb_raw.tt;
glmb_temp.w= zeros(length(cu),1);
glmb_temp.I= cell(length(cu),1);
glmb_temp.n= zeros(length(cu),1);
for hidx= 1:length(ic)
        glmb_temp.w(ic(hidx))= glmb_temp.w(ic(hidx))+glmb_raw.w(hidx);
        glmb_temp.I{ic(hidx)}= glmb_raw.I{hidx};
        glmb_temp.n(ic(hidx))= glmb_raw.n(hidx);
end
glmb_temp.cdn= glmb_raw.cdn;
end



function glmb_clean= clean_update(glmb_temp,filter)
%flag used tracks
usedindicator= zeros(length(glmb_temp.tt),1);
for hidx= 1:length(glmb_temp.w)
    usedindicator(glmb_temp.I{hidx})= usedindicator(glmb_temp.I{hidx})+1;
end
trackcount= sum(usedindicator>0);

%remove unused tracks and reindex existing hypotheses/components
newindices= zeros(length(glmb_temp.tt),1); newindices(usedindicator>0)= 1:trackcount;
glmb_clean.tt= glmb_temp.tt(usedindicator>0);
glmb_clean.w= glmb_temp.w;
for hidx= 1:length(glmb_temp.w)
    glmb_clean.I{hidx}= newindices(glmb_temp.I{hidx});
end
glmb_clean.n= glmb_temp.n;
glmb_clean.cdn= glmb_temp.cdn;

%resampling step for particle implementation
for tabidx= 1:length(glmb_clean.tt)
    neffsamptemp= 1/sum((glmb_clean.tt{tabidx}.w).^2);
    if neffsamptemp < filter.nth
        xtemptemp= glmb_clean.tt{tabidx}.x;
        wtemptemp= glmb_clean.tt{tabidx}.w;
        rspidx= randsample(length(wtemptemp),filter.npt,true,wtemptemp); %rspidx= resample(wtemptemp,filter.npt);
        glmb_clean.tt{tabidx}.x= xtemptemp(:,rspidx);
        glmb_clean.tt{tabidx}.w= ones(filter.npt,1)/filter.npt;
    end
end
end



function tt_lmb= glmb2lmb(glmb)

%find unique labels (with different possibly different association histories)
lmat= zeros(2,length(glmb.tt),1);
for tabidx= 1:length(glmb.tt)
    lmat(:,tabidx)= glmb.tt{tabidx}.l;
end
lmat= lmat';

[cu,~,ic]= unique(lmat,'rows'); cu= cu';

%initialize LMB struct
tt_lmb= cell(size(cu,2),1);
for tabidx=1:length(tt_lmb)
   tt_lmb{tabidx}.r= 0;
   tt_lmb{tabidx}.x= [];
   tt_lmb{tabidx}.w= [];
   tt_lmb{tabidx}.l= cu(:,tabidx);
end

%extract individual tracks
for hidx=1:length(glmb.w)
   for t= 1:glmb.n(hidx)
      trkidx= glmb.I{hidx}(t);
      newidx= ic(trkidx);
      tt_lmb{newidx}.x= cat(2,tt_lmb{newidx}.x,glmb.tt{trkidx}.x);
      tt_lmb{newidx}.w= cat(1,tt_lmb{newidx}.w,glmb.w(hidx)*glmb.tt{trkidx}.w);
   end
end

%extract existence probabilities and normalize track weights
for tabidx=1:length(tt_lmb)
   tt_lmb{tabidx}.r= sum(tt_lmb{tabidx}.w);
   tt_lmb{tabidx}.w= tt_lmb{tabidx}.w/tt_lmb{tabidx}.r;
end

end



function tt_lmb_out= clean_lmb(tt_lmb_in,filter)
%prune tracks with low existence probabilities
rvect= get_rvals(tt_lmb_in);
idxkeep= find(rvect > filter.track_threshold);
rvect= rvect(idxkeep);
tt_lmb_out= tt_lmb_in(idxkeep);

%enforce cap on maximum number of tracks
if length(tt_lmb_out) > filter.T_max
    [~,idxkeep]= sort(rvect,'descend');
    tt_lmb_out= tt_lmb_out(idxkeep);   
end

%cleanup tracks
for tabidx=1:length(tt_lmb_out)
    xtemptemp= tt_lmb_out{tabidx}.x;
    wtemptemp= tt_lmb_out{tabidx}.w;
    rspidx= randsample(length(wtemptemp),filter.npt,true,wtemptemp); %rspidx= resample(wtemptemp,filter.npt);
    tt_lmb_out{tabidx}.x= xtemptemp(:,rspidx);
    tt_lmb_out{tabidx}.w= ones(filter.npt,1)'/filter.npt;
end
end



function rvect= get_rvals(tt_lmb)                           %function to extract vector of existence probabilities from LMB track table
rvect= zeros(length(tt_lmb),1);
for tabidx=1:length(tt_lmb)
   rvect(tabidx)= tt_lmb{tabidx}.r; 
end
end



function [X,N,L]=extract_estimates(tt_lmb,model)
%extract estimates via MAP cardinality and corresponding tracks
rvect= get_rvals(tt_lmb); rvect= min(rvect,0.999); rvect= max(rvect,0.001);
cdn= prod(1-rvect)*esf(rvect./(1-rvect));
[~,mode] = max(cdn);
N = min(length(rvect),mode-1);
X= zeros(model.x_dim,N);
L= zeros(2,N);

[~,idxcmp]= sort(rvect,'descend');
for n=1:N
    [~,idxtrk]= max(tt_lmb{idxcmp(n)}.w);
    X(:,n)= tt_lmb{idxcmp(n)}.x*tt_lmb{idxcmp(n)}.w(:);
    L(:,n)= tt_lmb{idxcmp(n)}.l;
end
end



function display_diaginfo(tt_lmb,k,est,filter,T_predict,T_posterior,T_clean)
rvect= get_rvals(tt_lmb); rvect= min(rvect,0.999); rvect= max(rvect,0.001);
cdn= prod(1-rvect)*esf(rvect./(1-rvect));
eap= (0:(length(cdn)-1))*cdn(:);
var= (0:(length(cdn)-1)).^2*cdn(:)-((0:(length(cdn)-1))*cdn(:))^2;
if ~strcmp(filter.run_flag,'silence')
    disp([' time= ',num2str(k),...
        ' #eap cdn=' num2str(eap),...
        ' #var cdn=' num2str(var,4),...
        ' #est card=' num2str(est.N(k),4),...
        ' #trax pred=' num2str(T_predict,4),...
        ' #trax post=' num2str(T_posterior,4),...
        ' #trax updt=',num2str(T_clean,4)   ]);
end
end


function sorted_particles = kd_tree_matching(tt_lmb_update_x, tt_predict_x)
    % KD-Tree Based Nearest Neighbor Matching for Particle Sorting
    % Uses MATLAB's built-in knnsearch for fast nearest neighbor lookup.

    X1 = tt_lmb_update_x; % Updated particles
    X2 = tt_predict_x; % Predicted particles

    % Ensure feature dimensions match
    if size(X1, 1) ~= size(X2, 1)
        error('Dimension mismatch: Particles must have the same feature size.');
    end

    % Use KD-Tree nearest neighbor search
    assignment = knnsearch(X2', X1'); % Finds closest X2 point for each X1 point

    % Ensure unique assignments (optional)
    assignment = unique(assignment, 'stable');

    % Sort particles based on KD-Tree assignments
    sorted_particles = X1(:, assignment);

    % Ensure output type matches input
    sorted_particles = cast(sorted_particles, class(tt_lmb_update_x));
end


function sorted_particles = greedy_nearest_neighbor11(tt_lmb_update_x, tt_predict_x)
    % Greedy Nearest Neighbor Matching for Particle Sorting

    X1 = tt_lmb_update_x; % Updated particles
    X2 = tt_predict_x; % Predicted particles

    if size(X1, 1) ~= size(X2, 1)
        error('Dimension mismatch: Particles must have the same feature size.');
    end

    num_particles = size(X2, 2);  % Ensure the same number of particles

    % Compute pairwise distances
    cost_matrix = pdist2(X1', X2');

    % Initialize tracking variables
    assigned_particles = zeros(1, num_particles);  % Store selected indices
    available_indices = 1:size(X1, 2);  % Indices of available particles

    % Assign each predicted particle greedily
    for i = 1:num_particles
        % Find the closest available match for X2(:,i)
        [~, min_idx] = min(cost_matrix(:, i));

        % Assign this index
        assigned_particles(i) = available_indices(min_idx);

        % Remove this assigned index from availability to avoid duplicates
        available_indices(min_idx) = [];
        cost_matrix(min_idx, :) = Inf; % Prevent reassignment
    end

    % Sort particles based on assignments
    sorted_particles = X1(:, assigned_particles);

    % Ensure output size matches `tt_predict_x`
    if size(sorted_particles, 2) ~= num_particles
        error('Output size does not match tt_predict_x');
    end

    % Ensure output type matches input
    sorted_particles = cast(sorted_particles, class(tt_lmb_update_x));
end


function sorted_particles = sort_particles_hungarian(tt_lmb_update_x, tt_predict_x)
    % Hungarian Matching for Particle Sorting with Exact Size Matching

    % Extract matrices directly
    X1 = tt_lmb_update_x; % Updated particles
    X2 = tt_predict_x; % Predicted particles

    % Ensure they have compatible dimensions
    if size(X1, 1) ~= size(X2, 1)
        error('Dimension mismatch: Particles must have the same feature size.');
    end

    % Compute cost matrix (Euclidean distance)
    cost_matrix = pdist2(X1', X2'); % Correctly shaped cost matrix

    % Solve assignment using Hungarian algorithm with a large finite cost
    [row_idx, col_idx] = matchpairs(cost_matrix, 1e6); % Large but finite cost

    % Ensure row_idx contains only unique values
    row_idx = unique(row_idx, 'stable');

    % Restrict assignments to maintain correct size
    if length(row_idx) > size(X1, 2)
        row_idx = row_idx(1:size(X1, 2)); % Trim extra assignments
    elseif length(row_idx) < size(X1, 2)
        error('MatchPairs returned fewer assignments than expected.');
    end

    % Sort the updated particles based on matching order
    sorted_particles = X1(:, row_idx); 

    % Ensure output is of the same type
    sorted_particles = cast(sorted_particles, class(tt_lmb_update_x));
end



function sorted_particles = greedy_nearest_neighbor(particles_set1, particles_set2)
    % Function to sort particles from particles_set1 based on nearest neighbor matching with particles_set2
    
    num_particles1 = size(particles_set1, 2);
    num_particles2 = size(particles_set2, 2);
    
    if num_particles1 == 0 || num_particles2 == 0
        sorted_particles = [];
        return;
    end
    
    sorted_particles = zeros(size(particles_set1));
    
    % Track assigned particles
    assigned_indices = false(1, num_particles2);
    
    for i = 1:num_particles1
        min_dist = inf;
        best_match_idx = -1;
        
        % Find the nearest unassigned neighbor in particles_set2
        for j = 1:num_particles2
            if ~assigned_indices(j)
                dist = norm(particles_set1(:, i) - particles_set2(:, j));
                if dist < min_dist
                    min_dist = dist;
                    best_match_idx = j;
                end
            end
        end
        
        if best_match_idx ~= -1
            sorted_particles(:, i) = particles_set2(:, best_match_idx);
            assigned_indices(best_match_idx) = true;
        else
            sorted_particles(:, i) = particles_set1(:, i); % Keep the original if no match found
        end
    end
end
