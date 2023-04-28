clear
clc

run = 10;
npopArr = [40,80,100];                            % npopArr = [20,40];
for np = 1:length(npopArr)
    npop = npopArr(np);                                 % npop = 20,40
    
    maxitrArr = [10,30];                  % maxitrArr = [5,10];
    for mr = 1:length(maxitrArr)
        maxitr = maxitrArr(mr);                         % maxitr = 10,50,100,500
        
        nvarArr = [10,50,500,1000,1500,4000];        % nvarArr = [800,900];
        for nr = 1:length(nvarArr)
            nvar = nvarArr(nr);                         % nvar = 500,800,1000,1500,4000
            %%--------
            for maxrun = 1:run
                Elap = tic;
                time_allCode = tic;
                LB = -1+zeros(1,nvar);                 %lower bounds of variables
                UB = 1+zeros(1,nvar);                  %upper bounds of variables
                % Velocity Limits
                VelMax = 0.1*(UB-LB);
                VelMin = -VelMax;
                for m = 1:npop
                    for n = 1:nvar
                        particle(m).ArrayActThresh(1:nvar) = zeros(1,nvar); %#ok<SAGROW>
                        particle(m).Velocity(1:nvar) = zeros(1,nvar); %#ok<SAGROW>
                    end
                end
                finalmincost = zeros(1,maxitr);
                CostFunction=@(x, nvar) Exponential(x);        % Cost Function
                Func_time = tic;
                for m = 1:npop
                    for n = 1:nvar
                        particle(m).ArrayActThresh(1:nvar) = (LB(n)+rand(1,nvar)*(UB(n)-LB(n)));
                        particle(m).Velocity(1:nvar) = 0.1*particle(m).ArrayActThresh;
                    end
                end
                for m = 1:npop
                    % Evaluation
                    particle(m).cost=CostFunction(particle(m).ArrayActThresh);
                end
                time_func = toc(Func_time);
                particle_finally(1:npop) = particle(1:npop);
                for m = 1:npop
                    final_cost(m) = particle(m).cost; %#ok<SAGROW>
                end
                
                for x = 1:maxitr-1
                    [mincost, min_index] = min(final_cost);
                    Best = particle(min_index).ArrayActThresh;
                    [maxcost, max_index] = max(final_cost);
                    Worst = particle(max_index).ArrayActThresh;
                    finalmincost(1) = mincost;
                    pso_time = tic;
                    r = rand(npop,1);
                    r1 = rand(npop,1);
                    for m = 1:npop %loop start for updation of velocity and position i.e. ArraryActThersh
                        %% Constant inertia weight
                        w = 0.71;
                        c1 = 1.5;
                        c2 = 1.5;
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
                    end  % end of npop loop
                    time_pso = toc(pso_time);
                    
                    [min_final_cost, f_index] = min(final_cost);
                    finalmincost(x+1) = min_final_cost;
                    Best_particle = particle_finally(f_index).ArrayActThresh;
                end
                finalmincost = finalmincost';
                for i2 = 1:nvar
                    All(maxrun).all_Best_particle(i2) = Best_particle(i2);
                end
                All(maxrun).all_Best_particle_Cost = finalmincost;
                Elapsed(maxrun) = toc(Elap);%#ok<SAGROW>
                Serial_time_portionCode = time_func + time_pso;
                Serial_time_allCode = toc(time_allCode);
                
            end
            for e = 1:maxrun
                final_cost1(e) = mean(All(e).all_Best_particle_Cost); %#ok<SAGROW>
            end
            Result.avg_cost = mean(final_cost1);
            Result.min_cost = min(final_cost1);
            Result.max_cost = max(final_cost1);
            Result.deviation_cost = std(final_cost1);
            Result.avg_ET = mean(Elapsed);
            Result.Best_ET = min(Elapsed);
            Result.worst_ET = max(Elapsed);
        end
    end
end