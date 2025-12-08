function [eh,ev] = RefTarget(incidence,Param,facet_num)
% 根据材料选择发射率计算模型
% 输入：材料类型material，入射角incidence
% 输出：发射率的水平极化eh和垂直极化分量ev
switch Param.Target.material_type
    case 0
        [eh,ev] = RefCoherer(incidence,Param);
    case 1
        [eh,ev] = RefSingleConceal(incidence,Param);
    case 2
        [eh,ev] = RefTwoConceal(incidence,Param);
    case 3
        [eh,ev] = RefThreeConceal(incidence,Param);
    case 4
        [eh,ev] = RefFourConceal(incidence,Param);
    case 5
        [eh,ev] = RefDoubleTarget(incidence,Param,facet_num);
    otherwise
        error('This material is not supported');
end