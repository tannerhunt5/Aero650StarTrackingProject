// Fill out your copyright notice in the Description page of Project Settings.

using UnrealBuildTool;
using System.Collections.Generic;

public class StarTracking_Apr15Target : TargetRules
{
	public StarTracking_Apr15Target(TargetInfo Target) : base(Target)
	{
		Type = TargetType.Game;

		ExtraModuleNames.AddRange( new string[] { "StarTracking_Apr15" } );
	}
}
