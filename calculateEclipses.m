% Eclipse calculator script
% Predicts eclipses and syzygies (3+ body eclipses) between sets of planets 
% Takes a while to plot because it draws every X individually.

% (c) 2018 Sean Patrick 
% per MIT License
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function [finalPhases] = calculateEclipses (varargin)
close all
Epoch = 1000; %How long are we looking for eclipses? In this case, 100 days
% Units are whatever the orbit is in terms of. I.E. if orbital frequency is 1/19, an Epoch of 100 is 190 years.

% List of planets in order of distance from the sun separated by newlines
planetString = 'Abraxas \n Sseras \n Eldrick \n  Nex \n Brachine \n Folon \n Meket \n Zona \n Fornas \n Ariat \n';

% Orbital frequency in 1/(length of a year). Formatted on new lines per orbit ring.
ringFreqs = [1/19.2 1/54.6 1/11.2 1/15.0];

f = [ringFreqs(1) ringFreqs(1) ringFreqs(1)...
    ringFreqs(2) ringFreqs(2) ...
    ringFreqs(3) ringFreqs(3) ringFreqs(3) ...
    ringFreqs(4) ringFreqs(4)];

% Initial phase in radians. Planets with initial phase 0 are in eclipse at t=0.

if (~isempty(varargin))
    th = varargin{1};
else
    th = [0 pi/2 3*pi/4 0 pi  pi/4 2*pi/3 0 0 5*pi/8];
end


% The "shadow" of each planet, that is, how wide is the area it eclipses
% in terms of radian arc behind it? In order of distance from sun.

shadow = [0.3 0.2 0.1 0.2 0.4 0.2 0.3 0.1 0.1 0.1];

% The number of planets
Planets = length(f);

finalPhases = th;
for ind = 1:Planets
    finalPhases(ind) = orbitAngle(Epoch,f(ind),th(ind));
end
    
% Initializing...
eclipseMatrix = zeros(length(f),length(f),100);
eclipseCount = ones(length(f));
s=2;
syzygyMatrix = zeros(10,5);
lastSyzygy = [0 0 0];

% Main orbital loop
for t = 1:Epoch %At all t over the Epoch
    for j = 1:Planets % For each planet
        for k = 1:Planets %See if another planet lines up
       
            %Check if the current phase of two planets line up within their
            %shadow
            
            if (j<k) %Which planet is closer to the sun? ie. Whose shadow matters?
                jkShadow = shadow(j);
            else
                jkShadow = shadow(k);
            end
            
            if (abs(orbitAngle(t,f(j),th(j)) - orbitAngle(t,f(k),th(k))) < jkShadow && ...
               j ~= k && ...
               f(j) ~= f(k)) 
           
                % Store the time of the eclipse and increment the number of
                % eclipses that occurred.
                eclipseMatrix(j,k,eclipseCount(j,k)) = t;
                eclipseCount(j,k) = eclipseCount(j,k) + 1;
                
                
                %Detect syzygies
                for p = 1:Planets %Check if a 3rd planet is in line
                    for q = 1:Planets %Check if a 4th planet is in a line
                    
                        [planetsInvolved, planetOrder] = sort([j,k,p,q]);
                                  
                    
                        if (abs(orbitAngle(t,f(planetsInvolved(1)),th(planetsInvolved(1))) - orbitAngle(t,f(planetsInvolved(2)),th(planetsInvolved(2)))) < shadow(planetsInvolved(1)) && ...
                            abs(orbitAngle(t,f(planetsInvolved(2)),th(planetsInvolved(2))) - orbitAngle(t,f(planetsInvolved(3)),th(planetsInvolved(3)))) < shadow(planetsInvolved(2)) && ...
                            abs(orbitAngle(t,f(planetsInvolved(3)),th(planetsInvolved(3))) - orbitAngle(t,f(planetsInvolved(4)),th(planetsInvolved(4)))) < shadow(planetsInvolved(3)) && ...
                            length(unique([f(planetsInvolved(1)),f(planetsInvolved(2)),f(planetsInvolved(3)),f(planetsInvolved(4))])) == 4 && ...
                            ~isempty(setdiff(unique([j,k,p,q]),lastSyzygy))) %This is to make sure the same event doesn't get recorded twice per t

                            syzygyMatrix(s,:) = [t, j,k,p,q];
                            lastSyzygy = unique([j,k,p,q]); 
                            s = s+1;
                        end
                    
                    end %for q
                end %for p
            
            end %end syzygy check
                       
        end %end k loop
    end %end j loop
    

end %end t loop

eclipseCalendar = zeros(size(eclipseMatrix));

%Trim the zeros out
    for j = 1:length(f)
        for k = 1:length(f)
            eclipseDates = nonzeros(eclipseMatrix(j,k,:));
            if (~isempty(eclipseDates))
                eclipseCalendar(j,k,1:length(nonzeros(eclipseMatrix(j,k,:)))) = nonzeros(eclipseMatrix(j,k,:));
                eclipseCalendar(j,k,length(nonzeros(eclipseMatrix(j,k,:)))+1:end) = zeros(1,length(eclipseMatrix(j,k,:))-length(nonzeros(eclipseMatrix(j,k,:))));
            else
            end
        end
    end
eclipseCalendar(:,:,all(~eclipseCalendar,3)) = [];

%Plotting

%Color gradient used for the set of planets
gradient = 1/Planets:1/Planets:1;
colors = [gradient;zeros(1,length(gradient));1-gradient];

for p1 = 1:Planets
    figure()
    
    for p2 = 1:Planets
        
        %Don't waste time plotting within the same orbit. 
        if (f(p1) == f(p2))
            continue
        end
         
                %Trim zeros
                eclipseTimes = squeeze(eclipseCalendar(p1,p2,:));
                
                hold on
                
                %Mark an X at the other planets in an eclipse line
                h =plot(eclipseTimes, p2, 'x');  
                
                %Also mark an X at this planet
                h2 =plot(eclipseTimes, p1, 'x');  
                
                %Set the color of the X
                set(h,'Color',colors(:,p2))
                set(h2,'Color','k')
                
                axis([1 Epoch 0 Planets+1]);
                set(gca, 'XTick', Epoch/10:Epoch/10:Epoch)
                
                ytlblstr = sprintf(planetString);
                ytlbls = regexp(ytlblstr, '\n', 'split');
                ytlbls{end} = 'PLANETS';
                set(gca, 'YTick',1:length(ytlbls), 'YTickLabel',ytlbls)
               
                title(sprintf('Eclipse calendar for %s',ytlbls{p1}));
                xlabel('Hours') %I picked hours for my epoch unit, but it's arbitrary

        % Draw a vertical line where syzygies occurred
        for ind = 1:length(syzygyMatrix(:,1))
            if ( (syzygyMatrix(ind,2) == p1) || ...
                 (syzygyMatrix(ind,3) == p1) || ...
                 (syzygyMatrix(ind,4) == p1) || ...
                 (syzygyMatrix(ind,5) == p1) )
                
                line([syzygyMatrix(ind,1), syzygyMatrix(ind,1)],[0 Planets+1],'Color',[0.5 0.5 0.5]);   
                
                if (isempty(setdiff(unique([1,4,8,9]),unique(syzygyMatrix(ind,2:5)))))
                         %If it's the same set of planets as initial conditions
                    line([syzygyMatrix(ind,1), syzygyMatrix(ind,1)],[0 Planets+1],'Color',[1 0 0]); 
                end
            end
        end
        
    end
end

end

%Find the phase of a planet at time t given year length 1/f and
%initial phase th
function theta = orbitAngle(t,f,th)
theta = mod(2*pi*(f).*t + th,2*pi);
end