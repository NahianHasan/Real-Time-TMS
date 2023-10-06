function [Efields,Time,Memory,data_set_up_time,communication_time,GM_teid,GM_msh] = real_time_field_calculation_multiple_placements(Q,Ax,params,Transformations,memory_calculation)
    multiplying_factor = 6.6E7;
    Efields = zeros(size(Q,1),size(Transformations,3));
    g = params{31};
    NModes = params{29};
    num_simulations = size(Transformations,3);
    Time = []; Memory = [];
    if numel(g)>0
        overhead_1 = 0;
        for jx=1:num_simulations
            start_t = tic();
            Tr = gpuArray(single(Transformations(:,:,jx)));
            Tri = gpuArray(single(inv(Tr)));
            overhead_1 = overhead_1 + toc(start_t);
            [Efield_gpu,Tm,mem] = real_time_gpu(Q,Ax,params{1},params{2},params{3},params{4},params{5},...
                params{6},params{7},params{8},params{9},params{10},params{11},params{12},params{13},params{14},params{15},...
                params{16},params{17},params{18},params{19},params{20},params{21},params{22},params{23},params{24},...
                Tr,Tri,params{29},memory_calculation);
            Efields(:,jx) = gather(Efield_gpu);
            Time(jx,:,:) = [Tm(1),Tm(2)-Tm(1),Tm(3)-Tm(2),Tm(4)-Tm(3)];
            if memory_calculation
                Memory(jx) = mem;
            end
            disp(['GPU Running - Modes (w) = ',num2str(NModes),'--- ID=',num2str(jx),'--- Time=',num2str(Tm(4)*1000),' ms']);
        end
        overhead_1 = overhead_1/num_simulations;
        communication_time = params{33} + overhead_1;
        data_set_up_time = params{32} + overhead_1;
    else
        for jx=1:num_simulations
            [Efield_cpu,Tm] = real_time_cpu(Q,Ax,params{1},params{2},params{3},params{4},params{8},params{17},params{18},...
                                            params{19},params{20},params{21},params{22},Transformations(:,:,jx),params{29});
            disp(['CPU Running - Modes (w) = ',num2str(NModes),'--- ID=',num2str(jx),'--- Time=',num2str(Tm(4)*1000),' ms']);
            Efields(:,jx) = Efield_cpu;
            Time(jx,:,:) = [Tm(1),Tm(2)-Tm(1),Tm(3)-Tm(2),Tm(4)-Tm(3)];
            if memory_calculation
                T = whos;
                mem = 0;
                for ix=1:length(T)
                    mem = mem + T(ix).bytes/1024/1024;
                end
                Memory(jx) = mem;
            end
        end
    end
    GM_teid = params{27};
    GM_msh = params{28};
    Efields = Efields.*multiplying_factor;
end