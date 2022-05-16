function RecordMembrane = SetUpRecordMembrane(TimeVec,Membrane)

    RecordMembrane.Nodes = zeros(size(Membrane.Nodes,1),2,length(TimeVec),'single');
    RecordMembrane.Segments = Membrane.Segments;
    RecordMembrane.Springs = Membrane.Springs;
    RecordMembrane.TimeVec = TimeVec;   

end
