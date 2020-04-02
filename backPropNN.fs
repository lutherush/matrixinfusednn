let backProp (yVect : Vector<float>) (fwdPropEval : (NodesVectList * (Network))) : Network =
    let rec _calcPrevDelta (deltaLst : Vector<float> list) (layerLst: Network) (nodesInLayerLst: NodesVect list) : Vector<float> list =
        match layerLst with
        | [_] -> deltaLst
        | head :: head2 :: _ -> 
                    let interim = head.wMatx.RemoveRow(0) * deltaLst.Head
                    let theta' = nodesInLayerLst.Head.sVect.Map (Func<float, float> (getActFnDeriv (head2.actFn)))

                    //hamadard product
                    let prevDelta = interim.PointwiseMultiply(theta')
                    _calcPrevDelta (prevDelta :: deltaLst) (layerLst.Tail) (nodesInLayerLst.Tail)
        | [] -> [] //this will never hit if wMatx is checked initially
    
    let layerListUpdater (prevNodes : NodesVect) (currDelta : Vector<float>) (currLayer : Layer) : Layer =
        //must add 1.0 on top for correct shape.
        let dEdWMatx =  CreateMatrix.DenseOfColumnArrays( (Array.concat [ [|1.0|] ; prevNodes.xVect.ToArray() ])  ) * currDelta.ToRowMatrix()
        let newWMatx = currLayer.wMatx - (lr*dEdWMatx)
        {currLayer with wMatx=newWMatx}
        
    let xAndSLstRev, netRev = fwdPropEval
    let deltaL = lastLayerDeriv xAndSLstRev.Head.xVect xAndSLstRev.Head.sVect yVect netRev.Head.actFn
    let deltaRevLst = _calcPrevDelta [deltaL] netRev xAndSLstRev.Tail |> List.rev
    
    List.map3 layerListUpdater xAndSLstRev.Tail deltaRevLst netRev |> List.rev
