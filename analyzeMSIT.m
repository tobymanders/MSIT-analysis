function [behavioralStats,neuralStats] = analyzeMSIT(patientID,session,nevFile)
%UNTITLED2 Summary of this function goes here

[behavioralStats] = analyzeMSITbehavior(patientID,session,nevFile)

[neuralStats] = analyzeMSITunits (patientID,session,nevFile)


end

