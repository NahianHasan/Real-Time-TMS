function [Memory] = collect_memory(grid_spacing,Mode_Array,FEMORD,output_folders,coil_models)
    Memory_gpu = zeros([length(coil_models),length(Mode_Array),length(output_folders)]);
    Memory_cpu = zeros([length(coil_models),length(Mode_Array),length(output_folders)]);
    %collect memory
    for cx=1:length(coil_models)
        coil_model = coil_models{cx};
        for jx=1:length(output_folders)
            output_folder = output_folders{jx};
            for ix=1:length(Mode_Array)
                NModes = Mode_Array(ix);
                try
                    W = load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Memory_grid_',num2str(grid_spacing),'_Modes_',num2str(NModes),'_coil_',coil_model,'_gpu.mat']));
                    W1 = load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Memory_grid_',num2str(grid_spacing),'_Modes_',num2str(NModes),'_coil_',coil_model,'_cpu.mat']));
                catch
                    continue;
                end
                Memory_gpu(cx,ix,jx) = W.Memory/1024;%convert to GB
                Memory_cpu(cx,ix,jx) = W1.Memory/1024;%convert to GB
            end
        end
    end
    average_memory = mean(mean(Memory_gpu,3),1);
    max_memory = max(max(Memory_gpu,[],3),[],1);
    min_memory = min(min(Memory_gpu,[],3),[],1);
    Memory = struct();
    Memory.gpu1_memory_during_real_time = struct();
    Memory.gpu1_memory_during_real_time.memory = Memory_gpu;
    Memory.gpu1_memory_during_real_time.average_memory = average_memory;
    Memory.gpu1_memory_during_real_time.max_memory = max_memory;
    Memory.gpu1_memory_during_real_time.min_memory = min_memory;
    average_memory = mean(mean(Memory_cpu,3),1);
    max_memory = max(max(Memory_cpu,[],3),[],1);
    min_memory = min(min(Memory_cpu,[],3),[],1);
    Memory.cpu1_memory_during_real_time_in_cpu = struct();
    Memory.cpu1_memory_during_real_time_in_cpu.memory = Memory_cpu;
    Memory.cpu1_memory_during_real_time_in_cpu.average_memory = average_memory;
    Memory.cpu1_memory_during_real_time_in_cpu.max_memory = max_memory;
    Memory.cpu1_memory_during_real_time_in_cpu.min_memory = min_memory;
    Memory_cpu_gpu = Memory_cpu - Memory_gpu;
    average_memory = mean(mean(Memory_cpu_gpu,3),1);
    max_memory = max(max(Memory_cpu_gpu,[],3),[],1);
    min_memory = min(min(Memory_cpu_gpu,[],3),[],1);
    Memory.cpu1_memory_during_real_time_in_gpu = struct();
    Memory.cpu1_memory_during_real_time_in_gpu.memory = Memory_cpu_gpu;
    Memory.cpu1_memory_during_real_time_in_gpu.average_memory = average_memory;
    Memory.cpu1_memory_during_real_time_in_gpu.max_memory = max_memory;
    Memory.cpu1_memory_during_real_time_in_gpu.min_memory = min_memory;
end

