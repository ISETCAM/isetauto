function [ scene, mappings ] = MexximpRemodellerMultipleObj( scene, mappings, names, conditionValues, conditionNumber )

shadowDirection = eval(rtbGetNamedValue(names,conditionValues,'shadowDirection',[]));
cameraPosition = eval(rtbGetNamedValue(names,conditionValues,'cameraPosition',[]));
cameraLookAt = eval(rtbGetNamedValue(names,conditionValues,'cameraLookAt',[]));

cameraPTR = eval(rtbGetNamedValue(names,conditionValues,'cameraPTR','[0 0 0]'));

objMovementFile = rtbGetNamedValue(names,conditionValues,'objPosFile','');


%% Add a camera
scene = mexximpCentralizeCamera(scene);
lookUp = [0 0 -1];
cameraLookDir = cameraLookAt - cameraPosition;
cameraPTR = deg2rad(cameraPTR);

transformation = mexximpLookAt(1000*cameraPosition,1000*cameraLookAt,lookUp);
ptrTransform = mexximpPTR(cameraPTR(1), cameraPTR(2), cameraPTR(3), cameraLookDir, lookUp);


cameraId = strcmp({scene.rootNode.children.name},'Camera');
scene.rootNode.children(cameraId).transformation = transformation*mexximpTranslate(-1000*cameraPosition)*ptrTransform*mexximpTranslate(1000*cameraPosition);


%% Translate the objects

objects = loadjson(objMovementFile,'SimplifyCell',1);

for i=1:length(scene.rootNode.children)
    for o=1:length(objects)
        if isempty(strfind(scene.rootNode.children(i).name,objects(o).prefix)) == false
    
            position = eval(objects(o).position)*1000;
            orientation = eval(objects(o).orientation);
            
            scene.rootNode.children(i).transformation = scene.rootNode.children(i).transformation*...
           mexximpRotate([0 0 -1],deg2rad(orientation))*mexximpTranslate(position);
        end
       
   end
end




% Add directional light ('SunLight');
ambient = mexximpConstants('light');
ambient.position = [0 0 0]';
ambient.type = 'directional';
ambient.name = 'SunLight';
ambient.lookAtDirection = shadowDirection(:);
ambient.ambientColor = 10000*[1 1 1]';
ambient.diffuseColor = 10000*[1 1 1]';
ambient.specularColor = 10000*[1 1 1]';
ambient.constantAttenuation = 1;
ambient.linearAttenuation = 0;
ambient.quadraticAttenuation = 1;
ambient.innerConeAngle = 0;
ambient.outerConeAngle = 0;

scene.lights = [scene.lights, ambient];

ambientNode = mexximpConstants('node');
ambientNode.name = ambient.name;
ambientNode.transformation = eye(4);

scene.rootNode.children = [scene.rootNode.children, ambientNode];




end

