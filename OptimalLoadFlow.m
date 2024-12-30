% Run load flow with optimal photovoltaic DG size and placement
%global scenario

busdg=round(OV(1,1));
Pdg=OV(1,2);
Sdg=Pdg/PF;
%Sdg2=Pdg2/PF2;
%Sdg3=Pdg3/PF3;

Qdg=sqrt(Sdg^2-Pdg^2);
%Qdg2=sqrt(Sdg^2-Pdg2^2);
%Qdg3=sqrt(Sdg^2-Pdg3^2);
 

P(busdg,1)=P(busdg,1)-Pdg/(1000*MVAb);
Q(busdg,1)=Q(busdg,1)-Qdg/(1000*MVAb);
 
%P(busdg,1)=P(busdg,1)-Pdg2/(1000*MVAb);
%Q(busdg,1)=Q(busdg,1)-Qdg2/(1000*MVAb);

%P(busdg,1)=P(busdg,1)-Pdg3/(1000*MVAb);
%Q(busdg,1)=Q(busdg,1)-Qdg3/(1000*MVAb);

for s=1:10
for i=1:no
    nlc(i,1)=conj(complex(P(i,1),Q(i,1)))/(Vb(i,1)); % Find Ibr with root bus I=S*/V (At each branch)
end
nlc;
for i=1:br
    Ibr(i,1)=nlc(i+1,1); % Deleting root bus, and add to "Ibr"
end
Ibr;
xy=length(adjcb(1,:));
for i=br-1:-1:1
    for k=1:xy
        if adjcb(i,k)~=0
            u=adjcb(i,k);
            %Ibr(i,1)=nlc(i+1,1)+Ibr(k,1);
            Ibr(i,1)=Ibr(i,1)+Ibr(u,1); % Backward sweep, increment of I from end branch
        end
    end      
end
Ibr; % True value of I at respectively branch (After increment)
for i=2:no
      g=0;
      for a=1:b 
          if xy>1
            if adjcb(a,2)==i-1 
                u=adjcb(a,1);
                Vb(i,1)=((Vb(u,1))-((Ibr(i-1,1))*(complex((R(i-1,1)),X(i-1,1))))); % Backward sweep of voltage with first lateral connection
                g=1;
            end
            if adjcb(a,3)==i-1 
                u=adjcb(a,1);
                Vb(i,1)=((Vb(u,1))-((Ibr(i-1,1))*(complex((R(i-1,1)),X(i-1,1))))); % Backward sweep of voltage with second lateral connection
                g=1;
            end
          end
        end
        if g==0
            Vb(i,1)=((Vb(i-1,1))-((Ibr(i-1,1))*(complex((R(i-1,1)),X(i-1,1))))); % Backward sweep of voltage drop in normal case
        end
end
s=s+1;
end
nlc;
Ibr;
Vbp=[abs(Vb) angle(Vb)*180/pi];
 
for i=1:no
    Va(i,2:3)=Vbp(i,1:2);
end
for i=1:no
    Va(i,1)=i;
end
Va;
 
Ibrp=[abs(Ibr) angle(Ibr)*180/pi];
PL(1,1)=0;
QL(1,1)=0;
 
% Losses
for f=1:br
    Pl(f,1)=(Ibrp(f,1)^2)*R(f,1);
    Ql(f,1)=X(f,1)*(Ibrp(f,1)^2);
    PL(1,1)=PL(1,1)+Pl(f,1);
    QL(1,1)=QL(1,1)+Ql(f,1);
end

%Plosskw=(Pl)*100000;
%Qlosskw=(Ql)*100000;
%PL=(PL)*100000;
%QL=(QL)*100000;

Voltage = Vbp(:,1);
Angle = Vbp(:,2)*(pi/180);

