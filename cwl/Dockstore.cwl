#!/usr/bin/env cwl-runner
cwlVersion: v1.2
$namespaces:
      cwltool: http://commonwl.org/cwltool#
$graph:
  - class: Workflow
    id: main
    label: Workflow that executes the SBG/Sister Fractional Cover workflow
    
    requirements:
      InlineJavascriptRequirement: {}
      NetworkAccess:
        networkAccess: true
      StepInputExpressionRequirement: {}
      SubworkflowFeatureRequirement: {}

    inputs:
      parameters:
        type:
          fields:
            crid:
            - 'null'
            - string
            output_collection:
            - 'null'
            - string
            cores:
            - 'null'
            - int
            refl_scale:
            - 'null'
            -int
            normalization:
            - 'null'
            - string
            experimental:
            - 'null'
            - boolean
          type: record
      stage_in:
        type:
          fields:
            download_type: string
            downloading_keys:
            - string
            - 'null'
            edl_password:
            - string
            - 'null'
            edl_password_type:
            - string
            - 'null'
            edl_username:
            - string
            - 'null'
            stac_json: File
          type: record
      stage_out:
        type:
          fields:
            aws_access_key_id:
            - string
            - 'null'
            aws_region:
            - string
            - 'null'
            aws_secret_access_key:
            - string
            - 'null'
            aws_session_token:
            - string
            - 'null'
            collection_id:
            - string
            - 'null'
            staging_bucket:
            - string
            - 'null'
          type: record
      stage_aux_in:
        type:
          fields:
            download_type: string
            downloading_keys:
            - string
            - 'null'
            edl_password:
            - string
            - 'null'
            edl_password_type:
            - string
            - 'null'
            edl_username:
            - string
            - 'null'
            stac_json: File
          type: record
    outputs:
      stage_out_results:
        outputSource: stage_out/stage_out_results
        type: File
      stage_out_success:
        outputSource: stage_out/successful_features
        type: File
      stage_out_failures:
        outputSource: stage_out/failed_features
        type: File

    steps:
      process:
          run: "#unity-fractional-cover"
          in:
            #catalog_name: input_catalog_name
            crid:
              source: parameters
              valueFrom: $(self.crid)
            download_dir: stage_in/stage_in_download_dir
            aux_dir: stage_aux_in/stage_in_download_dir
            output_collection: 
              source: parameters
              valueFrom: $(self.output_collection)
            cores: 
              source: parameters
              valueFrom: $(self.cores)
            refl_scale: 
              source: parameters
              valueFrom: $(self.refl_scale)
            normalization: 
              source: parameters
              valueFrom: $(self.normalization)
            experimental: 
              source: parameters
              valueFrom: $(self.experimental)
          out: [process_output_dir]
      stage_in:
        in:
          download_type:
            source: stage_in
            valueFrom: $(self.download_type)
          downloading_keys:
            source: stage_in
            valueFrom: $(self.downloading_keys)
          edl_password:
            source: stage_in
            valueFrom: $(self.edl_password)
          edl_password_type:
            source: stage_in
            valueFrom: $(self.edl_password_type)
          edl_username:
            source: stage_in
            valueFrom: $(self.edl_username)
          stac_json:
            source: stage_in
            valueFrom: $(self.stac_json)
        out:
        - stage_in_collection_file
        - stage_in_download_dir
        run: stage_in.cwl
      stage_aux_in:
        in:
          download_type:
            source: stage_in
            valueFrom: $(self.download_type)
          downloading_keys:
            source: stage_in
            valueFrom: $(self.downloading_keys)
          edl_password:
            source: stage_in
            valueFrom: $(self.edl_password)
          edl_password_type:
            source: stage_in
            valueFrom: $(self.edl_password_type)
          edl_username:
            source: stage_in
            valueFrom: $(self.edl_username)
          stac_json:
            source: stage_aux_in
            valueFrom: $(self.stac_json)
        out:
        - stage_in_collection_file
        - stage_in_download_dir
        run: stage_in.cwl
      stage_out:
        in:
          aws_access_key_id:
            source: stage_out
            valueFrom: $(self.aws_access_key_id)
          aws_secret_access_key:
            source: stage_out
            valueFrom: $(self.aws_secret_access_key)
          aws_session_token:
            source: stage_out
            valueFrom: $(self.aws_session_token)
          collection_id:
            source: stage_out
            valueFrom: $(self.collection_id)
          output_dir: process/process_output_dir
          staging_bucket:
            source: stage_out
            valueFrom: $(self.staging_bucket)
        out:
        - stage_out_results
        - successful_features
        - failed_features
        run: stage_out.cwl
    
  - class: CommandLineTool
    id: "sbg-unity-fractional-cover"
    label: "sbg fractional cover"
    cwlVersion: v1.2
    
    requirements:
      - class: DockerRequirement
        dockerPull: "leebrian0731/sbg-unity-fractional-cover:"
    
    
    arguments:
    - $(inputs.download_dir.path)/$(inputs.catalog_name)
    - $(runtime.outdir)
    - $(inputs.cores)
    - $(inputs.refl_scale)
    - $(inputs.normalization)
    - $(inputs.crid)
    - $(inputs.experimental)
    - $(inputs.output_collection)
    inputs:
      catalog_name:
        type: string
        default: stage-in-results.json
      crid:
        default: '001'
        type: string
      download_dir:
        type:  Directory
      output_collection_name:
        default: 'SBG-L2-Fractional-Cover'
        type: string
      cores:
        default: 1
        type: int
      refl_scale:
        default: 0
        type: int
      normalization:
        default: 'none'
        type: string
      experimental:
        default: True
        type: boolean
      

    
    outputs:
      process_output_dir:
        outputBinding:
          glob: $(runtime.outdir)
        type: Directory
    
    baseCommand: ["python", "process.py"]
