/*
 * Copyright (C) 1997-2001 Id Software, Inc.
 * Copyright (C) 2011 Knightmare
 * Copyright (C) 2011 Yamagi Burmeister
 * Copyright (c) ZeniMax Media Inc.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 *
 * =======================================================================
 *
 * Prototypes for every mmove_t in the game.so.
 *
 * =======================================================================
 */

extern mmove_t tank_move_death ;
extern mmove_t tank_move_attack_chain ;
extern mmove_t tank_move_attack_post_rocket ;
extern mmove_t tank_move_attack_fire_rocket ;
extern mmove_t tank_move_attack_pre_rocket ;
extern mmove_t tank_move_attack_strike ;
extern mmove_t tank_move_attack_post_blast ;
extern mmove_t tank_move_reattack_blast ;
extern mmove_t tank_move_attack_blast ;
extern mmove_t tank_move_pain3 ;
extern mmove_t tank_move_pain2 ;
extern mmove_t tank_move_pain1 ;
extern mmove_t tank_move_stop_run ;
extern mmove_t tank_move_run ;
extern mmove_t tank_move_start_run ;
extern mmove_t tank_move_stop_walk ;
extern mmove_t tank_move_walk ;
extern mmove_t tank_move_start_walk ;
extern mmove_t tank_move_stand ;
extern mmove_t supertank_move_end_attack1 ;
extern mmove_t supertank_move_attack1 ;
extern mmove_t supertank_move_attack2 ;
extern mmove_t supertank_move_attack3 ;
extern mmove_t supertank_move_attack4 ;
extern mmove_t supertank_move_backward ;
extern mmove_t supertank_move_death ;
extern mmove_t supertank_move_pain1 ;
extern mmove_t supertank_move_pain2 ;
extern mmove_t supertank_move_pain3 ;
extern mmove_t supertank_move_turn_left ;
extern mmove_t supertank_move_turn_right ;
extern mmove_t supertank_move_forward ;
extern mmove_t supertank_move_run ;
extern mmove_t supertank_move_stand ;
extern mmove_t soldierh_move_death6 ;
extern mmove_t soldierh_move_death5 ;
extern mmove_t soldierh_move_death4 ;
extern mmove_t soldierh_move_death3 ;
extern mmove_t soldierh_move_death2 ;
extern mmove_t soldierh_move_death1 ;
extern mmove_t soldierh_move_duck ;
extern mmove_t soldierh_move_attack6 ;
extern mmove_t soldierh_move_attack4 ;
extern mmove_t soldierh_move_attack3 ;
extern mmove_t soldierh_move_attack2 ;
extern mmove_t soldierh_move_attack1 ;
extern mmove_t soldierh_move_pain4 ;
extern mmove_t soldierh_move_pain3 ;
extern mmove_t soldierh_move_pain2 ;
extern mmove_t soldierh_move_pain1 ;
extern mmove_t soldierh_move_run ;
extern mmove_t soldierh_move_start_run ;
extern mmove_t soldierh_move_walk2 ;
extern mmove_t soldierh_move_walk1 ;
extern mmove_t soldierh_move_stand3 ;
extern mmove_t soldierh_move_stand1 ;
extern mmove_t soldier_move_death6 ;
extern mmove_t soldier_move_death5 ;
extern mmove_t soldier_move_death4 ;
extern mmove_t soldier_move_death3 ;
extern mmove_t soldier_move_death2 ;
extern mmove_t soldier_move_death1 ;
extern mmove_t soldier_move_duck ;
extern mmove_t soldier_move_attack6 ;
extern mmove_t soldier_move_attack4 ;
extern mmove_t soldier_move_attack3 ;
extern mmove_t soldier_move_attack2 ;
extern mmove_t soldier_move_attack1 ;
extern mmove_t soldier_move_pain4 ;
extern mmove_t soldier_move_pain3 ;
extern mmove_t soldier_move_pain2 ;
extern mmove_t soldier_move_pain1 ;
extern mmove_t soldier_move_run ;
extern mmove_t soldier_move_start_run ;
extern mmove_t soldier_move_walk2 ;
extern mmove_t soldier_move_walk1 ;
extern mmove_t soldier_move_stand3 ;
extern mmove_t soldier_move_stand1 ;
extern mmove_t parasite_move_death ;
extern mmove_t parasite_move_break ;
extern mmove_t parasite_move_drain ;
extern mmove_t parasite_move_pain1 ;
extern mmove_t parasite_move_stop_walk ;
extern mmove_t parasite_move_start_walk ;
extern mmove_t parasite_move_walk ;
extern mmove_t parasite_move_stop_run ;
extern mmove_t parasite_move_start_run ;
extern mmove_t parasite_move_run ;
extern mmove_t parasite_move_stand ;
extern mmove_t parasite_move_end_fidget ;
extern mmove_t parasite_move_fidget ;
extern mmove_t parasite_move_start_fidget ;
extern mmove_t mutant_move_death2 ;
extern mmove_t mutant_move_death1 ;
extern mmove_t mutant_move_pain3 ;
extern mmove_t mutant_move_pain2 ;
extern mmove_t mutant_move_pain1 ;
extern mmove_t mutant_move_jump ;
extern mmove_t mutant_move_attack ;
extern mmove_t mutant_move_run ;
extern mmove_t mutant_move_start_walk ;
extern mmove_t mutant_move_walk ;
extern mmove_t mutant_move_idle ;
extern mmove_t mutant_move_stand ;
extern mmove_t medic_move_attackCable ;
extern mmove_t medic_move_attackBlaster ;
extern mmove_t medic_move_attackHyperBlaster ;
extern mmove_t medic_move_duck ;
extern mmove_t medic_move_death ;
extern mmove_t medic_move_pain2 ;
extern mmove_t medic_move_pain1 ;
extern mmove_t medic_move_run ;
extern mmove_t medic_move_walk ;
extern mmove_t medic_move_stand ;
extern mmove_t insane_move_struggle_cross ;
extern mmove_t insane_move_cross ;
extern mmove_t insane_move_crawl_death ;
extern mmove_t insane_move_crawl_pain ;
extern mmove_t insane_move_runcrawl ;
extern mmove_t insane_move_crawl ;
extern mmove_t insane_move_stand_death ;
extern mmove_t insane_move_stand_pain ;
extern mmove_t insane_move_run_insane ;
extern mmove_t insane_move_walk_insane ;
extern mmove_t insane_move_run_normal ;
extern mmove_t insane_move_walk_normal ;
extern mmove_t insane_move_down ;
extern mmove_t insane_move_jumpdown ;
extern mmove_t insane_move_downtoup ;
extern mmove_t insane_move_uptodown ;
extern mmove_t insane_move_stand_insane ;
extern mmove_t insane_move_stand_normal ;
extern mmove_t infantry_move_attack2 ;
extern mmove_t infantry_move_attack1 ;
extern mmove_t infantry_move_duck ;
extern mmove_t infantry_move_death3 ;
extern mmove_t infantry_move_death2 ;
extern mmove_t infantry_move_death1 ;
extern mmove_t infantry_move_pain2 ;
extern mmove_t infantry_move_pain1 ;
extern mmove_t infantry_move_run ;
extern mmove_t infantry_move_walk ;
extern mmove_t infantry_move_fidget ;
extern mmove_t infantry_move_stand ;
extern mmove_t hover_move_end_attack ;
extern mmove_t hover_move_attack1 ;
extern mmove_t hover_move_start_attack ;
extern mmove_t hover_move_backward ;
extern mmove_t hover_move_death1 ;
extern mmove_t hover_move_run ;
extern mmove_t hover_move_walk ;
extern mmove_t hover_move_forward ;
extern mmove_t hover_move_land ;
extern mmove_t hover_move_pain1 ;
extern mmove_t hover_move_pain2 ;
extern mmove_t hover_move_pain3 ;
extern mmove_t hover_move_takeoff ;
extern mmove_t hover_move_stop2 ;
extern mmove_t hover_move_stop1 ;
extern mmove_t hover_move_stand ;
extern mmove_t gunner_move_attack_grenade ;
extern mmove_t gunner_move_endfire_chain ;
extern mmove_t gunner_move_fire_chain ;
extern mmove_t gunner_move_attack_chain ;
extern mmove_t gunner_move_duck ;
extern mmove_t gunner_move_death ;
extern mmove_t gunner_move_pain1 ;
extern mmove_t gunner_move_pain2 ;
extern mmove_t gunner_move_pain3 ;
extern mmove_t gunner_move_runandshoot ;
extern mmove_t gunner_move_run ;
extern mmove_t gunner_move_walk ;
extern mmove_t gunner_move_stand ;
extern mmove_t gunner_move_fidget ;
extern mmove_t gladiator_move_death ;
extern mmove_t gladiator_move_pain_air ;
extern mmove_t gladiator_move_pain ;
extern mmove_t gladiator_move_attack_gun ;
extern mmove_t gladiator_move_attack_melee ;
extern mmove_t gladiator_move_run ;
extern mmove_t gladiator_move_walk ;
extern mmove_t gladiator_move_stand ;
extern mmove_t gladb_move_death ;
extern mmove_t gladb_move_pain_air ;
extern mmove_t gladb_move_pain ;
extern mmove_t gladb_move_attack_gun ;
extern mmove_t gladb_move_attack_melee ;
extern mmove_t gladb_move_run ;
extern mmove_t gladb_move_walk ;
extern mmove_t gladb_move_stand ;
extern mmove_t gekk_move_rduck ;
extern mmove_t gekk_move_lduck ;
extern mmove_t gekk_move_wdeath ;
extern mmove_t gekk_move_death4 ;
extern mmove_t gekk_move_death3 ;
extern mmove_t gekk_move_death1 ;
extern mmove_t gekk_move_pain2 ;
extern mmove_t gekk_move_pain1 ;
extern mmove_t gekk_move_pain ;
extern mmove_t gekk_move_attack ;
extern mmove_t gekk_move_leapatk2 ;
extern mmove_t gekk_move_leapatk ;
extern mmove_t gekk_move_attack2 ;
extern mmove_t gekk_move_attack1 ;
extern mmove_t gekk_move_spit ;
extern mmove_t gekk_move_run_start ;
extern mmove_t gekk_move_run ;
extern mmove_t gekk_move_walk ;
extern mmove_t gekk_move_idle2 ;
extern mmove_t gekk_move_idle ;
extern mmove_t gekk_move_swim_start ;
extern mmove_t gekk_move_swim_loop ;
extern mmove_t gekk_move_standunderwater ;
extern mmove_t gekk_move_stand ;
extern mmove_t flyer_move_loop_melee ;
extern mmove_t flyer_move_end_melee ;
extern mmove_t flyer_move_start_melee ;
extern mmove_t flyer_move_attack2 ;
extern mmove_t flyer_move_bankleft ;
extern mmove_t flyer_move_bankright ;
extern mmove_t flyer_move_defense ;
extern mmove_t flyer_move_pain1 ;
extern mmove_t flyer_move_pain2 ;
extern mmove_t flyer_move_pain3 ;
extern mmove_t flyer_move_rollleft ;
extern mmove_t flyer_move_rollright ;
extern mmove_t flyer_move_stop ;
extern mmove_t flyer_move_start ;
extern mmove_t flyer_move_run ;
extern mmove_t flyer_move_walk ;
extern mmove_t flyer_move_stand ;
extern mmove_t floater_move_run ;
extern mmove_t floater_move_walk ;
extern mmove_t floater_move_pain3 ;
extern mmove_t floater_move_pain2 ;
extern mmove_t floater_move_pain1 ;
extern mmove_t floater_move_death ;
extern mmove_t floater_move_attack3 ;
extern mmove_t floater_move_attack2 ;
extern mmove_t floater_move_attack1 ;
extern mmove_t floater_move_activate ;
extern mmove_t floater_move_stand2 ;
extern mmove_t floater_move_stand1 ;
extern mmove_t flipper_move_death ;
extern mmove_t flipper_move_attack ;
extern mmove_t flipper_move_pain1 ;
extern mmove_t flipper_move_pain2 ;
extern mmove_t flipper_move_start_run ;
extern mmove_t flipper_move_walk ;
extern mmove_t flipper_move_run_start ;
extern mmove_t flipper_move_run_loop ;
extern mmove_t flipper_move_stand ;
extern mmove_t fixbot_move_weld_end ;
extern mmove_t fixbot_move_weld ;
extern mmove_t fixbot_move_weld_start ;
extern mmove_t fixbot_move_attack2 ;
extern mmove_t fixbot_move_laserattack ;
extern mmove_t fixbot_move_attack1 ;
extern mmove_t fixbot_move_start_attack ;
extern mmove_t fixbot_move_backward ;
extern mmove_t fixbot_move_death1 ;
extern mmove_t fixbot_move_run ;
extern mmove_t fixbot_move_walk ;
extern mmove_t fixbot_move_forward ;
extern mmove_t fixbot_move_land ;
extern mmove_t fixbot_move_pain3 ;
extern mmove_t fixbot_move_painb ;
extern mmove_t fixbot_move_paina ;
extern mmove_t fixbot_move_takeoff ;
extern mmove_t fixbot_move_turn ;
extern mmove_t fixbot_move_roamgoal ;
extern mmove_t fixbot_move_pickup ;
extern mmove_t fixbot_move_stand2 ;
extern mmove_t fixbot_move_stand ;
extern mmove_t fixbot_move_landing ;
extern mmove_t chick_move_start_slash ;
extern mmove_t chick_move_end_slash ;
extern mmove_t chick_move_slash ;
extern mmove_t chick_move_end_attack1 ;
extern mmove_t chick_move_attack1 ;
extern mmove_t chick_move_start_attack1 ;
extern mmove_t chick_move_duck ;
extern mmove_t chick_move_death1 ;
extern mmove_t chick_move_death2 ;
extern mmove_t chick_move_pain3 ;
extern mmove_t chick_move_pain2 ;
extern mmove_t chick_move_pain1 ;
extern mmove_t chick_move_walk ;
extern mmove_t chick_move_run ;
extern mmove_t chick_move_start_run ;
extern mmove_t chick_move_stand ;
extern mmove_t chick_move_fidget ;
extern mmove_t brain_move_run ;
extern mmove_t brain_move_attack4 ;
extern mmove_t brain_move_attack3 ;
extern mmove_t brain_move_attack2 ;
extern mmove_t brain_move_attack1 ;
extern mmove_t brain_move_death1 ;
extern mmove_t brain_move_death2 ;
extern mmove_t brain_move_duck ;
extern mmove_t brain_move_pain1 ;
extern mmove_t brain_move_pain2 ;
extern mmove_t brain_move_pain3 ;
extern mmove_t brain_move_defense ;
extern mmove_t brain_move_walk1 ;
extern mmove_t brain_move_idle ;
extern mmove_t brain_move_stand ;
extern mmove_t boss5_move_end_attack1 ;
extern mmove_t boss5_move_attack1 ;
extern mmove_t boss5_move_attack2 ;
extern mmove_t boss5_move_attack3 ;
extern mmove_t boss5_move_attack4 ;
extern mmove_t boss5_move_backward ;
extern mmove_t boss5_move_death ;
extern mmove_t boss5_move_pain1 ;
extern mmove_t boss5_move_pain2 ;
extern mmove_t boss5_move_pain3 ;
extern mmove_t boss5_move_turn_left ;
extern mmove_t boss5_move_turn_right ;
extern mmove_t boss5_move_forward ;
extern mmove_t boss5_move_run ;
extern mmove_t boss5_move_stand ;
extern mmove_t makron_move_attack5 ;
extern mmove_t makron_move_attack4 ;
extern mmove_t makron_move_attack3 ;
extern mmove_t makron_move_sight ;
extern mmove_t makron_move_death3 ;
extern mmove_t makron_move_death2 ;
extern mmove_t makron_move_pain4 ;
extern mmove_t makron_move_pain5 ;
extern mmove_t makron_move_pain6 ;
extern mmove_t makron_move_walk ;
extern mmove_t makron_move_run ;
extern mmove_t makron_move_stand ;
extern mmove_t jorg_move_end_attack1 ;
extern mmove_t jorg_move_attack1 ;
extern mmove_t jorg_move_start_attack1 ;
extern mmove_t jorg_move_attack2 ;
extern mmove_t jorg_move_death ;
extern mmove_t jorg_move_pain1 ;
extern mmove_t jorg_move_pain2 ;
extern mmove_t jorg_move_pain3 ;
extern mmove_t jorg_move_end_walk ;
extern mmove_t jorg_move_walk ;
extern mmove_t jorg_move_start_walk ;
extern mmove_t jorg_move_run ;
extern mmove_t jorg_move_stand ;
extern mmove_t boss2_move_death ;
extern mmove_t boss2_move_pain_light ;
extern mmove_t boss2_move_pain_heavy ;
extern mmove_t boss2_move_attack_rocket ;
extern mmove_t boss2_move_attack_post_mg ;
extern mmove_t boss2_move_attack_mg ;
extern mmove_t boss2_move_attack_pre_mg ;
extern mmove_t boss2_move_run ;
extern mmove_t boss2_move_walk ;
extern mmove_t boss2_move_fidget ;
extern mmove_t boss2_move_stand ;
extern mmove_t berserk_move_death2 ;
extern mmove_t berserk_move_death1 ;
extern mmove_t berserk_move_pain2 ;
extern mmove_t berserk_move_pain1 ;
extern mmove_t berserk_move_attack_strike ;
extern mmove_t berserk_move_attack_club ;
extern mmove_t berserk_move_attack_spike ;
extern mmove_t berserk_move_run1 ;
extern mmove_t berserk_move_walk ;
extern mmove_t berserk_move_stand_fidget ;
extern mmove_t berserk_move_stand ;
