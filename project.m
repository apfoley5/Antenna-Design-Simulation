clear;
clc;

%arrays initialization
A = zeros(181); %holds theta values
B = zeros(181); %holds abs(f(theta)) values
C = zeros(181,360); %holds abs(g(theta,psi) values
omega = zeros(181,360);

%initial values
m = 2; %2
n = 3; %3
psi1 = 30; %30
psi2 = 258; %258

%alpha, d and N values
alpha = (psi1+psi2)/2
d = (psi2-psi1)/(4*pi)
N = m*(n-1) + 1
%filling matrix B
for i = 1:181 %theta
    A(i,:) = i-1;
    B(i,:) = abs(arrayFactor(psi1, psi2, m, n, (i-1)));
    C_yz(i,:) = abs(elementPattern(i-1,90));
    C_xz(i,:) = abs(elementPattern(i-1,0));
    
    for j = 1:360 %phi
        C(i,j) = abs(elementPattern(i-1,j-1));
    end
    totalPattern_xz(i,:) = B(i)*C_xz(i);
    totalPattern_yz(i,:) = B(i)*C_yz(i);
    totalPattern = B*C;
    totalPattern = totalPattern(:,2);
end

%finding the normalized total pattern
max_xz = 0;
max_yz = 0;
max = 0;
for value1 = 1:181
    if max_xz < totalPattern_xz(value1)
        max_xz = totalPattern_xz(value1);
    end
    if max_yz < totalPattern_yz(value1)
        max_yz = totalPattern_yz(value1);
    end
    if max < totalPattern(value1)
        max = totalPattern(value1);
    end
end

%finding the normalized total pattern
normalizedTotalPattern_xz = 20*log10(totalPattern_xz/max_xz);
normalizedTotalPattern_yz = 20*log10(totalPattern_yz/max_yz);
normalizedTotalPattern = 20*log10(totalPattern/max);

%fixing the normalized total pattern
for value2 = 1:181
    if normalizedTotalPattern_xz(value2) <= -50
        normalizedTotalPattern_xz(value2) = -50;
    end
    if normalizedTotalPattern_yz(value2) <= -50
        normalizedTotalPattern_yz(value2) = -50;
    end
end

figure(1)
polarplot(deg2rad(A),normalizedTotalPattern_xz)
title('xz plane');
rlim([-50 0]);
figure(2)
polarplot(deg2rad(A),normalizedTotalPattern_yz)
title('yz plane');
rlim([-50 0]);

%calculating Omega_a
const = (pi/180)^2;
omegaSum = 0;
for i4 = 1:181
    for j2 = 1:360
        x = B(i4);
        y = C(i4,j2);
        omega(i4,j2) = (const)*(((x*y)^2)*sind(i4-1));
    end
end
%this is the sum of all of the omega values
omegaSum = nansum(nansum(omega));
%directivity value
Directivity = (4*pi)/omegaSum;
Directivity_dB = 10*log10(Directivity)

%Half-Power Beamwidth for yz-plane
%first find proper index I1
halfpower_var_yz = (1/(sqrt(2))*B(181));
halfpower_var_2_yz = 100;
index_hp_yz = 0;
for p = 1:181
    if abs(halfpower_var_yz - B(p)) <= halfpower_var_2_yz
        halfpower_var_2_yz = abs(halfpower_var_yz - B(p));
        index_hp_yz = p;
    end
end

%half power value for yz-plane
HP_yz = 2*(180-A(index_hp_yz))

%Half-Power Beamwidth for xz-plane
halfpower_var_xz = (1/(sqrt(2))*B(181)*C(181,1));
halfpower_var_2_xz = 100;
index_hp_xz = 0;
for p2 = 1:181
    if abs(halfpower_var_xz - (B(p2)*C(p2,1))) <= halfpower_var_2_xz
        halfpower_var_2_xz = abs(halfpower_var_xz - (B(p2)*C(p2,1)));
        index_hp_xz = p2;
    end
end

%half power value for xz-plane
HP_xz = 2*(180-A(index_hp_xz))

%side lobe level in xz-plane
%finding the second largest lobe maximum in the normalized total pattern
secondLargestValue = 0;
for value1_2 = 2:180
    if( (normalizedTotalPattern_xz(value1_2) > normalizedTotalPattern_xz(value1_2-1)) & (normalizedTotalPattern_xz(value1_2) > normalizedTotalPattern_xz(value1_2+1)) )
        secondLargestValue = normalizedTotalPattern_xz(value1_2);
    end
end

SLL = secondLargestValue
