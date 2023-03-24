#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

    Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

        // TODO Project 3-2: Part 1
        // Implement MirrorBSDF
        *pdf = 1;
        reflect(wo, wi);
        return reflectance / abs_cos_theta(*wi);
    }

    void MirrorBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Mirror BSDF"))
        {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            ImGui::TreePop();
        }
    }

// Microfacet BSDF //

    double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
        return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
    }

    double MicrofacetBSDF::D(const Vector3D h) {
        // TODO Project 3-2: Part 2
        // Compute Beckmann normal distribution function (NDF) here.
        // You will need the roughness alpha.
        return exp(-pow(tan(acos(h.z)), 2) / pow(alpha, 2)) / (PI * pow(alpha, 2) * pow(h.z, 4));
    }

    Vector3D MicrofacetBSDF::F(const Vector3D wi) {
        // TODO Project 3-2: Part 2
        // Compute Fresnel term for reflection on dielectric-conductor interface.
        // You will need both eta and etaK, both of which are Vector3D.
        Vector3D tmp1 = eta * eta + k * k;
        Vector3D tmp2 = 2 * eta * abs_cos_theta(wi);
        double tmp3 = pow(abs_cos_theta(wi), 2);
        Vector3D Rs = (tmp1 - tmp2 + tmp3) / (tmp1 + tmp2 + tmp3);
        Vector3D Rp = (tmp1 * tmp3 - tmp2 + 1) / (tmp1 * tmp3 + tmp2 + 1);
        return (Rs + Rp) / 2;
    }

    Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
        // TODO Project 3-2: Part 2
        // Implement microfacet model here.
        if (wi.z > 0 && wo.z > 0) {
            Vector3D h = (wi + wo) / (wi + wo).norm();
            return F(wi) * G(wo, wi) * D(h) / 4 / wo.z / wi.z;
        }
        return Vector3D();
    }

    Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
        // TODO Project 3-2: Part 2
        // *Importance* sample Beckmann normal distribution function (NDF) here.
        // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
        //       and return the sampled BRDF value.

        // *wi = cosineHemisphereSampler.get_sample(pdf);
        // return MicrofacetBSDF::f(wo, *wi);

        Vector2D r = sampler.get_sample();
        double r1 = r[0], r2 = r[1];

        double theta_h = atan(sqrt(-1 * pow(alpha, 2) * log(1 - r1)));
        double phi_h = 2 * PI * r2;

        Vector3D h = Vector3D(sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h));
        *wi = -wo + 2 * dot(wo, h) * h;
        wi->normalize();
        if (wi->z < 0) {
            *pdf = 0;
            return Vector3D();
        }

        double p_theta = (2 * sin(theta_h)) / (pow(alpha, 2) * pow(cos(theta_h), 3)) * exp(-1 * pow(tan(theta_h) / alpha, 2));
        double p_phi = 1 / (2 * PI);
        double pw = p_theta * p_phi / sin(theta_h);
        *pdf = pw / (4 * dot(*wi, h));

        return MicrofacetBSDF::f(wo, *wi);
    }

    void MicrofacetBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Micofacet BSDF"))
        {
            DragDouble3("eta", &eta[0], 0.005);
            DragDouble3("K", &k[0], 0.005);
            DragDouble("alpha", &alpha, 0.005);
            ImGui::TreePop();
        }
    }

// Refraction BSDF //

    Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
        // TODO Project 3-2: Part 1
        // Implement RefractionBSDF
        if (!refract(wo, wi, ior)) {
            return Vector3D();
        }
        *pdf = 1;
        double eta = (wo.z > 0) ? 1. / ior : ior;
        return transmittance / abs_cos_theta(*wi) / pow(eta, 2);
    }

    void RefractionBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Refraction BSDF"))
        {
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

// Glass BSDF //

    Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

        // TODO Project 3-2: Part 1
        // Compute Fresnel coefficient and either reflect or refract based on it.

        // compute Fresnel coefficient and use it as the probability of reflection
        // - Fundamentals of Computer Graphics page 305
        double eta = (wo.z > 0) ? 1.0 / ior : ior;
        double totalInternalReflection = 1 - pow(eta, 2) * (1 - pow(wo.z, 2));
        // only reflect
        if (totalInternalReflection < 0) {
            *pdf = 1;
            reflect(wo, wi);
            return reflectance / abs_cos_theta(*wi);
        }

        // both reflect and refract
        double R0 = pow((ior - 1) / (ior + 1), 2);
        double R = R0 + (1 - R0) * pow(1 - abs(wo.z), 5);
        if (coin_flip(R)) {
            *pdf = R;
            reflect(wo, wi);
            return R * reflectance / abs_cos_theta(*wi);
        }
        else {
            *pdf = 1 - R;
            refract(wo, wi, ior);
            return (1 - R) * transmittance / abs_cos_theta(*wi) / pow(eta, 2);
        }
        return Vector3D();
    }

    void GlassBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Refraction BSDF"))
        {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

    void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

        // TODO Project 3-2: Part 1
        // Implement reflection of wo about normal (0,0,1) and store result in wi.
        // (xi, yi, zi) = (-xo, -yo, zo)
        wi->x = -wo.x;
        wi->y = -wo.y;
        wi->z = wo.z;
    }

    bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {

        // TODO Project 3-2: Part 1
        // Use Snell's Law to refract wo surface and store result ray in wi.
        // Return false if refraction does not occur due to total internal reflection
        // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
        // ray entering the surface through vacuum.
        double eta = (wo.z > 0) ? 1.0 / ior : ior;
        double totalInternalReflection = 1 - pow(eta, 2) * (1 - pow(wo.z, 2));
        if (totalInternalReflection < 0) {
            return false;
        }
        wi->x = -eta * wo.x;
        wi->y = -eta * wo.y;
        wi->z = (wo.z > 0) ? -sqrt(totalInternalReflection) : sqrt(totalInternalReflection);
        wi->normalize();
        return true;

    }

} // namespace CGL
