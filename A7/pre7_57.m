function Z = pre7_57(model,S,steps)
    X = Polyhedron('lb',model.x.min, 'ub',model.x.max);
    U = Polyhedron('lb',model.u.min, 'ub',model.u.max);
    Z = S;
    for i=1:steps
       R = model.reachableSet('X',Z, 'U',U, 'N',1, 'direction','backward');
       Z = X.intersect(R);
    end
end