%% Main program
%% ///////////////////////////////////////////////////////////////////
% Comments:
% 	Collaborator: Taku Nonomura, Yuji Saito, Keigo Yamada, 
%                 Kumi Nakai, Takayuki Nagata
% 	Last modified: 2021/01/05
% Nomenclature:
% - Scalars
%   n : Number of degrees of freedoaaam of spatial POD modes (state dimension)
%   p : Number of sensors
%   r : Number of rank for truncated POD
%   m : Number of snaphot (temporal dimension)
%   s : Number of randomized space (sketch size)
% - Matrices
% 	X : Supervising data matrix
% 	Y : Observation matrix
% 	H : Sparse sensor location matrix
% 	U : Spatial POD modes matrix
% 	C : Measurement matrix
% 	Z : POD mode amplitude matrix
%% ===================================================================

clear; close all;
warning('off','all')

%% Selection of Problems ============================================
% num_problem=1; % //Randomized sensor problem//
  num_problem=2; % //NOAA-SST//
% !<NOAA-SST> It takes a long time to obtain the solution in the convex 
% !<NOAA-SST> approximation method and the convex method is commented out 
% !<NOAA_SST> as default setting for reduction of demo time.
%
%% Parameters =======================================================
r = 10;
% pmin = 1;
% pinc = 1;
% pmax = 1;
% ps   = pmin:pinc:pmax;
ps = 1;%[5 8 10];
num_ave = 3;%200; % Number of iteration for averaging operation
CNT = 0; % Counter
maxiteration = 200; % Max iteration for convex approximation
% //Randomized sensor problem//
n = 100;%2000;
s = floor(n/10);
% //NOAA-SST//
m = 52*10; % 10-years (52weeks/year)
num_video = 2; % maxmum: m

%% Preparation of output directories ================================
workdir   = ('../work');
videodir  = [workdir,'/video'];
sensordir = [workdir,'/sensor_location'];
mkdir(workdir);
mkdir(videodir);
mkdir(sensordir);

%% Randomized sensor problem ========================================
if num_problem == 1

    %% Sensor selection =============================================
    for p = ps
        CNT = CNT+1;
        text = [ num2str(p),' sensor selection started --->' ];
        disp(text);

        %% Average loop =============================================
        for w=1:1:num_ave

            %% Preprocess for Randomized problem ====================
            U = randn(n,r);

            %% Random selection -------------------------------------
            [time_rand(CNT,w+1), H_rand, sensors_rand] = F_sensor_random(n,p);
            det_rand (CNT,w+1) = F_calc_det  (p,H_rand,U);
            tr_rand  (CNT,w+1) = F_calc_trace(p,H_rand,U);
            eig_rand (CNT,w+1) = F_calc_eigen(p,H_rand,U);

            %% D-optimality - Convex---------------------------------
            %!! This is very time consuming proceduce, We do not recommend to try this
            % [time_DC(CNT,w+1), H_DC, sensors_DC, ...
            %  NT_TOL_cal_DC(CNT,w+1), iter_DC(CNT,w+1)] ...
            %  = F_sensor_DC(U,p,maxiteration);
            % det_DC (CNT,w+1) = F_calc_det  (p,H_DC,U);
            % tr_DC  (CNT,w+1) = F_calc_trace(p,H_DC,U);
            % eig_DC (CNT,w+1) = F_calc_eigen(p,H_DC,U);
            %!! I recommend you use the following dummy values
            %   if you do not need the solution in the convex approximation in NOAA-SST.        
            time_DC(CNT,w+1) = time_rand(CNT,w+1);
            det_DC (CNT,w+1) = det_rand (CNT,w+1);
            tr_DC  (CNT,w+1) = tr_rand  (CNT,w+1);
            eig_DC (CNT,w+1) = eig_rand (CNT,w+1);
            H_DC=H_rand;
            sensors_DC=sensors_rand;
            NT_TOL_cal_DC(CNT,w+1)=0;
            iter_DC(CNT,w+1)=0;
            
            %% Randomized Subspace Newton Convex Method (RSNC) -------
            [time_RSNC(CNT,w+1), H_RSNC, sensors_RSNC, ~] ...
            = F_sensor_RSNC(U,p,maxiteration,s);
            det_RSNC (CNT,w+1) = F_calc_det  (p,H_RSNC,U);
            tr_RSNC  (CNT,w+1) = F_calc_trace(p,H_RSNC,U);
            eig_RSNC (CNT,w+1) = F_calc_eigen(p,H_RSNC,U);
            
            %% Customized Randomized Subspace Newton Convex Method (CRSNC)
            [time_CRSNC(CNT,w+1), H_CRSNC, sensors_CRSNC, ~] ...
            = F_sensor_CRSNC(U,p,maxiteration,s);
            det_CRSNC (CNT,w+1) = F_calc_det  (p,H_RSNC,U);
            tr_CRSNC  (CNT,w+1) = F_calc_trace(p,H_RSNC,U);
            eig_CRSNC (CNT,w+1) = F_calc_eigen(p,H_RSNC,U);

        end
        
        %% Averaging ================================================
        [ time_rand, det_rand, tr_rand, eig_rand ]...
        = F_data_ave1( CNT, num_ave, time_rand, det_rand, tr_rand, eig_rand );
        [ time_DC, det_DC, tr_DC, eig_DC ]...
        = F_data_ave1( CNT, num_ave, time_DC, det_DC, tr_DC, eig_DC );
        [ time_RSNC, det_RSNC, tr_RSNC, eig_RSNC ]...
        = F_data_ave1( CNT, num_ave, time_RSNC, det_RSNC, tr_RSNC, eig_RSNC );
        [ time_RSNC, det_CRSNC, tr_CRSNC, eig_CRSNC ]...
        = F_data_ave1( CNT, num_ave, time_CRSNC, det_CRSNC, tr_CRSNC, eig_CRSNC );
        NT_TOL_cal_DC(CNT,1)=mean(NT_TOL_cal_DC(CNT,2:w+1));
        iter_DC(CNT,1)=mean(iter_DC(CNT,2:w+1));
        
        %% Sensor location ==========================================
        sensor_memo = zeros(p,4);
        sensor_memo(1:p,1) = sensors_rand(1:p)';
        sensor_memo(1:p,2) = sensors_DC(1:p);
        sensor_memo(1:p,3) = sensors_RSNC(1:p);
        sensor_memo(1:p,4) = sensors_CRSNC(1:p);
        filename = [workdir, '/sensors_p_', num2str(p), '.mat'];
        save(filename,'sensor_memo');

        text = [ '---> ', num2str(p), ' sensor selection finished!' ];
        disp(text);
    end
end

%% NOAA-SST =========================================================
if num_problem == 2

    %% Preprocces for NOAA-SST ======================================
    text='Readinng/Arranging a NOAA-SST dataset';
    disp(text);
    [Lat, Lon, time, mask, sst]...
    = F_pre_read_NOAA_SST( ['sst.wkmean.1990-present.nc'], ['lsmask.nc'] );
    [Uorg, Sorg, Vorg, Xorg, meansst, n] = F_pre_SVD_NOAA_SST(m, time, mask, sst);
    F_map_original(num_video, Xorg, meansst, mask, time, videodir);
    [U, Error_ave_pod, Error_std_pod]...
    = F_pre_truncatedSVD(r, Xorg, Uorg, Sorg, Vorg, num_video, meansst, mask, time, m, videodir);
    Error_ave_pod = repmat( Error_ave_pod , size(ps,2) );
    text='Complete Reading/Arranging a NOAA-SST dataset!';
    disp(text);
    s = floor(n/10);

    %% Sensor selection =============================================
    for p = ps
        CNT = CNT+1;
        text = [ num2str(p),' sensor selection started --->' ];
        disp(text);

        %% Random selection -----------------------------------------
        % Average loop
        for w=1:1:num_ave
            [time_rand(CNT,w+1), H_rand, sensors_rand] = F_sensor_random(n,p);
            det_rand(CNT,w+1) = F_calc_det  (p,H_rand,U);
            tr_rand (CNT,w+1) = F_calc_trace(p,H_rand,U);
            eig_rand(CNT,w+1) = F_calc_eigen(p,H_rand,U);
            [Zestimate_rand, Error_rand(CNT,w+1), Error_std_rand(CNT,w+1)] ...
            = F_calc_error(m, Xorg, U, H_rand);
        end
        % Averaging
        [ time_rand, det_rand, tr_rand, eig_rand, Error_rand, Error_std_rand ]...
        = F_data_ave2( CNT, num_ave, time_rand, det_rand, tr_rand, eig_rand, Error_rand, Error_std_rand );

        %% D-optimality - Convex-------------------------------------
        %!! This is very time consuming proceduce, We do not recommend to try this
        % [time_DC(CNT,1), H_DC, sensors_DC, ...
        %  NT_TOL_cal_DC(CNT,1), iter_DC(CNT,1)] ...
        %  = F_sensor_DC(U,p,maxiteration);
        % det_DC(CNT,1) = F_calc_det(p,H_DC,U);
        % tr_DC(CNT,1)  = F_calc_trace(p,H_DC,U);
        % eig_DC(CNT,1) = F_calc_eigen(p,H_DC,U);
        %!! I recommend you use the following dummy values 
        %   if you do not need the solution in the convex approximation in NOAA-SST.
        time_DC(CNT,1) = time_rand(CNT,1);
        det_DC (CNT,1) = det_rand (CNT,1);
        tr_DC  (CNT,1) = tr_rand  (CNT,1);
        eig_DC (CNT,1) = eig_rand (CNT,1);
        H_DC=H_rand;
        sensors_DC=sensors_rand;
        NT_TOL_cal_DC(CNT,w+1)=0;
        iter_DC(CNT,w+1)=0;
        %!!
        [Zestimate_DC, Error_DC(CNT,1), Error_std_DC(CNT,1)] ...
        = F_calc_error(m, Xorg, U, H_DC);
        NT_TOL_cal_DC(CNT,1)=mean(NT_TOL_cal_DC(CNT,2:w+1));
        iter_DC(CNT,1)=mean(iter_DC(CNT,2:w+1));
        
        %% Randomized Subspace Newton Convex Method (RSNC) -------
        [time_RSNC(CNT,1), H_RSNC, sensors_RSNC, ~] ...
        = F_sensor_RSNC(U,p,maxiteration,s);
        det_RSNC (CNT,1) = F_calc_det  (p,H_RSNC,U);
        tr_RSNC  (CNT,1) = F_calc_trace(p,H_RSNC,U);
        eig_RSNC (CNT,1) = F_calc_eigen(p,H_RSNC,U);
        [Zestimate_RSNC, Error_RSNC(CNT,1), Error_std_RSNC(CNT,1)] ...
        = F_calc_error(m, Xorg, U, H_RSNC);
        
        %% Customized Randomized Subspace Newton Convex Method (CRSNC)
        [time_CRSNC(CNT,1), H_CRSNC, sensors_CRSNC, ~] ...
        = F_sensor_CRSNC(U,p,maxiteration,s);
        det_CRSNC (CNT,1) = F_calc_det  (p,H_RSNC,U);
        tr_CRSNC  (CNT,1) = F_calc_trace(p,H_RSNC,U);
        eig_CRSNC (CNT,1) = F_calc_eigen(p,H_RSNC,U);
        [Zestimate_CRSNC, Error_CRSNC(CNT,1), Error_std_CRSNC(CNT,1)] ...
        = F_calc_error(m, Xorg, U, H_CRSNC);

        %% Sensor location ==========================================
        sensor_memo = zeros(p,4);
        sensor_memo(1:p,1) = sensors_rand(1:p)';
        sensor_memo(1:p,2) = sensors_DC(1:p);
        sensor_memo(1:p,3) = sensors_RSNC(1:p);
        sensor_memo(1:p,4) = sensors_CRSNC(1:p);
        filename = [workdir, '/sensors_p_', num2str(p), '.mat'];
        save(filename,'sensor_memo');

        %% Video ====================================================
        name='rand';
        F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
                        sensors_rand, Zestimate_rand, name, videodir, sensordir)
        name='DC';
        F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
                        sensors_DC, Zestimate_DC, name, videodir, sensordir)
        name='RSNC';
        F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
                        sensors_RSNC, Zestimate_RSNC, name, videodir, sensordir)
        name='CRSNC';
        F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
                        sensors_CRSNC, Zestimate_CRSNC, name, videodir, sensordir)
        
        text = [ '---> ', num2str(p), ' sensor selection finished!' ];
        disp(text);
    end
end

%% Data organization ================================================
% Arrange
[time_all] = F_data_arrange1( ps,   CNT, time_rand, time_DC,...
                              time_RSNC, time_CRSNC );
[det_all]  = F_data_arrange1( ps,   CNT, det_rand,  det_DC,...
                              det_RSNC,  det_CRSNC );
[tr_all]   = F_data_arrange1( ps,   CNT, tr_rand,   tr_DC,...
                              tr_RSNC,   tr_CRSNC );
[eig_all]  = F_data_arrange1( ps,   CNT, eig_rand,  eig_DC,...
                              eig_RSNC,  eig_CRSNC );
if num_problem == 2
    [Error] = F_data_arrange2( ps,          CNT, ...
                               Error_rand,  Error_std_rand, ...
                               Error_DC,    Error_std_DC,   ...
                               Error_RSNC,  Error_std_RSNC, ...
                               Error_CRSNC, Error_std_CRSNC,...
                               Error_ave_pod );
end
[log_DC] = F_data_arrange3( ps, CNT, NT_TOL_cal_DC, iter_DC );

% Normalize
[Normalized_det] = F_data_normalize( ps, CNT,  det_rand, det_DC,...
                                     det_RSNC, det_CRSNC );
[Normalized_tr]  = F_data_normalize( ps, CNT,  tr_rand,  tr_DC,...
                                     tr_RSNC,  tr_CRSNC );
[Normalized_eig] = F_data_normalize( ps, CNT,  eig_rand, eig_DC,...
                                     eig_RSNC, eig_CRSNC );

%% Save =============================================================
cd(workdir)
save('time.mat','time_all');
save('det.mat','det_all');
save('trace.mat','tr_all');
save('eigen.mat','eig_all');
save('Normalized_det.mat','Normalized_det');
save('Normalized_trace.mat','Normalized_tr');
save('Normalized_eigen.mat','Normalized_eig');
save('time_rand.mat','time_rand');
save('det_rand.mat','det_rand');
save('trace_rand.mat','tr_rand');
save('eigen_rand.mat','eig_rand');
if num_problem == 2
    save('Error.mat','Error');
    save('Error_rand.mat','Error_rand');
end
save('log_DC.mat','log_DC');

warning('on','all')
disp('Congratulations!');
cd ../src
%% ///////////////////////////////////////////////////////////////////
%% Main program end