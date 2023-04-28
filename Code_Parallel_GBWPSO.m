clear
clc
%% population parameter
run = 10;
npopArr = [40,80,100];                  % npopArr = [40,100,200,500];
for np = 1:length(npopArr)
    npop = npopArr(np);             % npop = 40,100,200,500
    
    maxitrArr = [10,30];             % maxitrArr = [10,50,100,500];
    for mr = 1:length(maxitrArr)
        maxitr = maxitrArr(mr);     % maxitr = 10,50,100,500
        
        nvarArr = [10,50,500,1000,1500,4000];   % nvarArr = [500,800,1000,1500,4000];
        for nr = 1:length(nvarArr)
            nvar = nvarArr(nr);     % nvar = 500,800,1000,1500,4000
            
            final_cost1 = zeros(1, run);
            ET = zeros(1, run);
            All = struct('all_Gbest', zeros(1, nvar),'all_Gbest_Cost', zeros(1, maxitr));
            All_time = tic;     % for time*****************************
            for maxrun = 1:run
                CostFunction = @(x, nVar) Exponential(x);        % Cost Function
                empty_particle = struct('ArrayActThresh', zeros(1, nvar),'Velocity', zeros(1, nvar));
                particle = repmat(empty_particle,npop,1);
                LB = -1+zeros(1,nvar); %lower bounds of variables
                UB = 1+zeros(1,nvar); %upper bounds of variables
                % Velocity Limits
                VelMax = 0.1*(UB-LB);
                VelMin = -VelMax;
                finalmincost = zeros(1,maxitr);
                Par_func = tic;     % for time*****************************
                parfor m = 1:npop
                    particle(m).ArrayActThresh = LB+rand(1,nvar).*(UB-LB);
                    particle(m).Velocity = 0.1*particle(m).ArrayActThresh;
                    % Evaluation
                    particle(m).cost = CostFunction(particle(m).ArrayActThresh);
                    final_cost(m) = particle(m).cost;
                end
                %                 Func_Par_time = toc(Par_func);     % for time*****************************
                %                 Best(1:npop) = particle(1:npop);  %  old Best
                for x = 1:maxitr-1
                    particle_finally(1:npop) = particle(1:npop);    %  particle_finally is that particle which will update finally in the algo. process
                    [mincost, min_index] = min(final_cost);
                    Best = particle(min_index).ArrayActThresh;
                    [maxcost, max_index] = max(final_cost);
                    Worst = particle(max_index).ArrayActThresh;
                    finalmincost(1) = mincost;
                    c1 = 1.5;
                    c2 = 1.5;
                    %% Constant inertia weight
                    w = 0.71;
                    r = rand(npop,1);
                    r1 = rand(npop,1);
                    parfor m = 1:npop %loop start for updation of velocity and position i.e. ArraryActThersh
                        % Update velocity
                        particle(m).Velocity = w * particle(m).Velocity + c1* r(m) * (Best - particle(m).ArrayActThresh) + c2 * r1(m) * (Worst - particle(m).ArrayActThresh);
                        % Apply Velocity Limits
                        particle(m).Velocity = max(particle(m).Velocity,VelMin);
                        particle(m).Velocity = min(particle(m).Velocity,VelMax);
                        % Update position
                        particle(m).ArrayActThresh = particle(m).ArrayActThresh + particle(m).Velocity;
                        % Velocity Mirror Effect
                        IsOutside = (particle(m).ArrayActThresh < LB | particle(m).ArrayActThresh > UB);
                        particle(m).Velocity(IsOutside) = -particle(m).Velocity(IsOutside);
                        % Apply Position Limits
                        particle(m).ArrayActThresh = max(particle(m).ArrayActThresh,LB);
                        particle(m).ArrayActThresh = min(particle(m).ArrayActThresh,UB);
                        % evaluation of cost
                        particle(m).cost = CostFunction(particle(m).ArrayActThresh);
                        cost(m) = particle(m).cost;
                        %--------------------------------------------
                        if (cost(m) < final_cost(m))  %old Best
                            final_cost(m) = cost(m);
                            particle_finally(m) = particle(m);  %old Best
                        end
                        % % % % %                         p(x,m) = Par.toc;
                    end  % end of npop loop
                    [min_final_cost, f_index] = min(final_cost);
                    finalmincost(x+1) = min_final_cost;
                    Best_particle = particle_finally(f_index).ArrayActThresh;
                    %                     pso_Par_time = toc(P_time);     % for time*****************************
                    %         w = w*wdamp;
                end
                % % % % %                 par_time = p;
                finalmincost = finalmincost';
                for i2 = 1:nvar
                    All(maxrun).all_Best_particle(i2) = Best_particle(i2);
                end
                All(maxrun).all_Best_particle_Cost = finalmincost;
                ET(maxrun) = toc(All_time);     % for time*****************************
            end
            %             paralleltime = Func_Par_time + pso_Par_time;     % for time*****************************
            for e = 1:maxrun
                final_cost1(e) = mean(All(e).all_Best_particle_Cost);
            end
            Result.avg_cost = mean(final_cost1);
            Result.min_cost = min(final_cost1);
            Result.max_cost = max(final_cost1);
            Result.deviation_cost = std(final_cost1);
            Result.avg_ET = mean(ET);
            Result.Best_ET = min(ET);
            Result.worst_ET = max(ET);
            
        end
    end
end